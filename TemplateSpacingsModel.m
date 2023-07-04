% Takes fit coeffs for template spacings over i and Tobs and fits them over
% fdot

clear;
clc;
close all;

ParsevalSNR = 1;

%% Read in data
fitDegree = pi;

if (ParsevalSNR == 1)
    filename = sprintf('TemplateSpacingFitCoeffsfitDegree_%fP.csv', fitDegree);
else 
    filename = sprintf('TemplateSpacingFitCoeffsfitDegree_%f.csv', fitDegree);
end

Data = readmatrix(filename);
fdotvec_sig = Data(:, 1);
data = Data(:, 2:2:end);
uncertainties = Data(:, 3:2:end);
NParams = length(data(1, :));

%% Plot Data
for k = 1:NParams
    figure(k);
    sk = scatter(log10(abs(fdotvec_sig)), data(:, k), '*k');
    titleString = sprintf('Template Spacing Fit Parameter a_{%d} vs fdot', k);
    title(titleString);
    xlabel('Base 10 log of fdot (log(Hz/s))');
    ylabel('Template Spacing Fit Param');
    sk.LineWidth = 3;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    grid on;
end

%% Fit Data

model_func = @(a, x) a(1).*log10(abs(x)).^(a(2));
fitCoeffs = zeros(1, NParams);
dCoeffs = zeros(1, NParams);
fitExps = zeros(1, NParams);
dExps = zeros(1, NParams);
fdot_sigp = linspace(min(fdotvec_sig), max(fdotvec_sig), 1000);

for k = 1:NParams
    [fitParamfitParams, dParams] = fitChiSquare(log10(abs(fdotvec_sig)), data(:, k), model_func, -rand(1, 2), zeros(size(fdotvec_sig)), uncertainties(:, k));
    %{
    fitCoeffs(k) = fitParamfitParams(1);
    dCoeffs(k) = dParams(1).d;
    fitExps(k) = fitParamfitParams(2);
    dExps = dParams(2).d;
    %}

    figure(k);
    hold on;
    pk = plot(log10(abs(fdot_sigp)), model_func(fitParamfitParams, log10(abs(fdot_sigp))), '-r');
    legend('Calculated Params from const fdot', 'fit', 'Location', 'north');
    hold off;
end

