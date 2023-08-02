% Takes fit params from fitSNRgTemplateSpacings.m and fits them against
% fdot

clear;
clc;
close all;

ParsevalSNR = 1;
TS = 'f'; % Choose f or fdot

%% Read in template spacing fit coeffs
fitDegree = 1;

if (strcmp(TS, 'f'))
    if (ParsevalSNR == 1)
        filename = sprintf('fTemplateSpacingFitCoeffsfitDegree_%fP.csv', fitDegree);
    else 
        filename = sprintf('fTemplateSpacingFitCoeffsfitDegree_%f.csv', fitDegree);
    end
else
    if (ParsevalSNR == 1)
        filename = sprintf('fdotTemplateSpacingFitCoeffsfitDegree_%fP.csv', fitDegree);
    else 
        filename = sprintf('fdotTemplateSpacingFitCoeffsfitDegree_%f.csv', fitDegree);
    end
end

FitCoeffs1 = readmatrix(filename);
fdotvec_sig = FitCoeffs1(:, 1);
fitCoeffs1 = FitCoeffs1(:, 2:2:end);
uncertainties1 = FitCoeffs1(:, 3:2:end);
NParams = length(fitCoeffs1(1, :));

%% Plot template spacing fit coeffs
for k = 1:NParams
    figure;
    ek = errorbar(log10(abs(fdotvec_sig)), fitCoeffs1(:, k), uncertainties1(:, k), '.k');
    titleString = sprintf('Template Spacing Fit Parameter a_{%d} vs fdot', k);
    title(titleString);
    xlabel('Base 10 log of fdot (log(Hz/s))');
    ylabel('Template Spacing Fit Param');
    ek.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    grid on;
end

%% Fit template spacing fit coeffs

lin_log_model = @(a, x) a(1).*log10(abs(x)) + a(2);
fdot_sigp = linspace(min(fdotvec_sig), max(fdotvec_sig), 1000);
fitParamfitParams = cell(NParams);

for k = 1:NParams
    [fitParamfitParams{k}, dParamParams] = fitChiSquare(fdotvec_sig, fitCoeffs1(:, k), lin_log_model, -rand(1, 2), zeros(size(fdotvec_sig)), uncertainties1(:, k));

    figure(k);
    hold on;
    pk = plot(log10(abs(fdot_sigp)), lin_log_model(fitParamfitParams{k}, fdot_sigp), '-r');
    legend('Calculated Params from const fdot', 'fit', 'Location', 'north');
    hold off;
end
