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
missing_ratio = 0.2;
outlier_ratio = 0.2;
mag_times = 5;

out_mag = max(max(max(double(Y))))*mag_times;
Omega = make_pattern(Ysz, missing_ratio);
O = make_outlier(Ysz, outlier_ratio, out_mag) .* Omega;


%% Example. Proposed: SOFIA (Online)
clear opts;
opts.lambda1    = 0.001;
opts.lambda2    = 0.001;
opts.lambda3    = 10;
opts.mu         = 0.1;
opts.phi        = 0.01;
opts.maxEpoch   = 300;
opts.tol        = 1e-3;
    
results = sofia_wrapper(Y, O, Omega, R, m, cycles, 0, opts);


%% Experiments
%% 1. NRE (Normalized Residual Error)
interval = 15;

figure;
semilogy(results.time(1:interval:end), results.nre(1:interval:end), '-r', 'LineWidth', 1); hold on;

grid on;
ax = gca;
l = legend('SOFIA');
title(['missing ratio: ', num2str(missing_ratio), ', outlier ratio: ', num2str(outlier_ratio), ', outlier magnitude: ', num2str(out_mag)]);
xlabel('Data Stream Time Index','FontName','Arial','FontSize',12,'FontWeight','bold');
ylabel('Normalized Residual Error (NRE)','FontName','Arial','FontSize',12,'FontWeight','bold');

%% 2. ART (Average Running Time)
figure;
art = [results.art_update];
hb = bar(art);
set(gca,'YScale','log');
l = legend('SOFIA');
xlabel('Datasets','FontName','Arial','FontSize',13,'FontWeight','bold');
ylabel('Average Running Time (ART)','FontName','Arial','FontSize',13,'FontWeight','bold');

%% 3. RAE (Running Average Error)
figure;
rae = [results.rae];
hb = bar(rae);
l = legend('SOFIA');
xlabel('Datasets','FontName','Arial','FontSize',13,'FontWeight','bold');
ylabel('Running Average Error (RAE)','FontName','Arial','FontSize',13,'FontWeight','bold');
