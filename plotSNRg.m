%% Plot SNRg vs i for Multiple Trials and Various Observation Times

clear;
clc;
close all;

%% Read in Data 

fdot_sig = -5.e-5;

if (fdot_sig == -5.e-9)
    Tobsvec = [8., 12., 16., 20., 24., 28., 32., 36., 40.]; 
else
    Tobsvec = [4., 16./3., 8., 12., 16., 20., 24., 28., 32.]; 
end
Nsegvec = zeros(size(Tobsvec));
SNRg_i = cell(1, length(Tobsvec));
ivec = cell(1, length(Tobsvec));

for k = 1:length(Tobsvec)
    filename = sprintf('SNRgfdot_%0.es-2Tobs_%0.fhrP.csv', fdot_sig, Tobsvec(k));
    SNRg_i{k} = csvread(filename);
    ivec{k} = 1:length(SNRg_i{k}(1, :));
    Nsegvec(k) = ivec{k}(end);
end
Ntrials = length(SNRg_i{1}(:, 1));

%% Plot data

% Create legend strings
legendStrings = cell(1, Ntrials+1);
for l = 1:Ntrials
    legendStrings{l} = ['Trial ', num2str(l)];
end
legendStrings{11} = 'Average';

% Create figure
for k = 1:length(Tobsvec)
    figure(k)
    hold on;
    for l = 1:Ntrials
        skl = scatter(ivec{k}, SNRg_i{k}(l, :), '.');
        skl.MarkerEdgeColor = rand(1, 3);
    end
    pk = plot(ivec{k}, mean(SNRg_i{k}), '-k');
    pk.LineWidth = 3;
    title(['SNRg_i vs i for Tobs = ', num2str(Tobsvec(k)), ' hr']);
    xlabel('i');
    ylabel('SNRg');
    legend(legendStrings, 'Location', 'Northwest');
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    grid on; 
    hold off;
end

%% Fit Data for Each Tobs

% Create two fit models
powfunc = @(a, x) log(a(1)) + log(x).*a(2);

% Calculate fit parameters
powFitParams = cell(1, length(Tobsvec));
for k = 1:length(Tobsvec)
    powFitParams{k} = fitChiSquare(ivec{k}, mean(log(SNRg_i{k})), powfunc, [1, 0.5],  zeros(size(ivec{k})), 10.*std(log(SNRg_i{k})).^2);
end

%% Plot Data and Fits for Each Tobs

figure(length(Tobsvec)+1)
for k = 1:length(Tobsvec)
    subplot(3, 3, k)
    hold on;
    pk21 = plot(ivec{k}, mean(SNRg_i{k}), '-k');
    pk22 = plot(ivec{k}, exp(powfunc(powFitParams{k}, ivec{k})), '-r');
    pk23 = plot(ivec{k}, mean(SNRg_i{k}(:, end))/sqrt(Nsegvec(k)).*sqrt(ivec{k}), '--b');
    title(['SNRg_i vs i for Tobs = ', num2str(Tobsvec(k)), ' hr']);
    xlabel('i');
    ylabel('SNRg');
    l1 = legend('mean values', ['power fit (exponent of ', num2str(powFitParams{k}(2)), ')'],  'sqrt function through origin and last point', 'Location', 'Southeast');
    l1.FontSize = 11;
    pk21.LineWidth = 3;
    pk22.LineWidth = 2;
    pk23.LineWidth = 1;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    grid on; 
    hold off;
end

%% Plot All Data

figure(length(Tobsvec)+2)
hold on;
for k = 1:length(Tobsvec)
    pk31 = scatter3(log10(Tobsvec(k)).*ones(size(mean(SNRg_i{k}))), log10(ivec{k}), log10(mean(SNRg_i{k})), '*b');
end
hold off;
title('SNRg over i and Tobs (log-log-log)');
xlabel('log(Tobs) (log(hr))');
ylabel('log(i)');
zlabel('log(SNR)');
grid on;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 16;
