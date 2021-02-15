clc; clear;

addpath('data/');
addpath('sofia/');
addpath('tensor_toolbox-v3.1/');

rng('default');
rng(1);

%% load tensor
load('network.mat');
Y = tenfun(@log2, T+1);


%% set up parameters
Ysz = size(Y);
N = ndims(Y);
R = 5;
m = 168;
cycles = 3;
outlier_ratio = 0.2;
mag = 5;

out_mag = max(max(max(double(Y))))*mag;
O = make_outlier(Ysz, outlier_ratio, out_mag); 

forecast_steps = 200;

%% Example. Proposed: SOFIA 
clear opts;
opts.lambda1    = 0.001;
opts.lambda2    = 0.001;
opts.lambda3    = 10;
opts.mu         = 0.1;
opts.phi        = 0.01;
opts.maxEpoch   = 300;
opts.tol        = 1e-3;

%% 1. missing ratio = 0%
rng(1);
missing_ratio = 0;
Omega = make_pattern(Ysz, missing_ratio); 

results_0 = sofia_wrapper(Y, O, Omega, R, m, cycles, forecast_steps, opts);


%% 2. missing ratio = 30%
rng(1);
missing_ratio = 0.3;
Omega = make_pattern(Ysz, missing_ratio); 

results_30 = sofia_wrapper(Y, O, Omega, R, m, cycles, forecast_steps, opts);


%% 2. missing ratio = 50%
rng(1);
missing_ratio = 0.5;
Omega = make_pattern(Ysz, missing_ratio); 

results_50 = sofia_wrapper(Y, O, Omega, R, m, cycles, forecast_steps, opts);


%% 2. missing ratio = 70%
rng(1);
missing_ratio = 0.7;
Omega = make_pattern(Ysz, missing_ratio); 

results_70 = sofia_wrapper(Y, O, Omega, R, m, cycles, forecast_steps, opts);

%% Experiments
%% NRE (Normalized Residual Error)
figure;
plot(results_0.forecast_time, results_0.forecast_nre, '-r', 'LineWidth', 1); hold on;
plot(results_30.forecast_time, results_30.forecast_nre, '-m', 'LineWidth', 1); hold on;
plot(results_50.forecast_time, results_50.forecast_nre, '-c', 'LineWidth', 1); hold on;
plot(results_70.forecast_time, results_70.forecast_nre, '-k', 'LineWidth', 1); hold on;
l = legend('SOFIA missing=0%','SOFIA missing=30%','SOFIA missing=50%','SOFIA missing=70%');
title(['Forecasting Error, outlier ratio: ', num2str(outlier_ratio), ', outlier magnitude: ', num2str(out_mag)]);
xlabel('Data Stream Time Index','FontName','Arial','FontSize',12,'FontWeight','bold');
ylabel('Normalized Residual Error (NRE)','FontName','Arial','FontSize',12,'FontWeight','bold');


%% RAE or AFE (Running Average Error = Average Forecasting Error)
f = figure;
xsize = 4.5;
ysize = 4;

rae = [results_0.forecast_rae, results_30.forecast_rae, results_50.forecast_rae, results_70.forecast_rae];
hb = bar(rae);
hb.FaceColor = 'flat';

hb.CData(1,:) = [0.89, 0.1, 0.11];
hb.CData(2,:) = [1, 0.33, 0.33];
hb.CData(3,:) = [1, 0.55, 0.55];
hb.CData(4,:) = [1, 0.7, 0.7];

set(gca, 'FontSize', 15);
set(gca,'XTickLabel',{'0%','30%','50%','70%'});
xlab = xlabel('Missing Ratio','FontName','Helvetica','FontSize',20,'FontWeight','normal');
ylab = ylabel('Average Forecasting Error','FontName','Helvetica','FontSize',20,'FontWeight','normal');

set(gca, 'box', 'off');
set(gcf,'units','inches', 'position', [5,5,xsize,ysize]);

