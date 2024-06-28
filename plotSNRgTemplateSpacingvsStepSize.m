%% Plotting SNRg Template Spacing vs fdot Step Size used to calculate template spacing

clear;
clc;
close all;

ParsevalSNR = 1;
TS = 'fdot'; % Choose f or fdot

%% Read Data
% Choose frequency derivative (Hz/s), observation time (hr)
fdot_sig = -5.e-6;
Tobs_hr = 24.;

% Read in data and assign to variables
if (strcmp(TS, 'f'))
    if (ParsevalSNR == 1)
        filename = sprintf('SNRgfTemplateSpacingvsStepSizefdot_%0.eTobs_%0.fP.csv', fdot_sig, Tobs_hr);
    else
        filename = sprintf('SNRgfTemplateSpacingvsStepSizefdot_%0.eTobs_%0.f.csv', fdot_sig, Tobs_hr);
    end
else
    if (ParsevalSNR == 1)
        filename = sprintf('SNRgfdotTemplateSpacingvsStepSizefdot_%0.eTobs_%0.fP.csv', fdot_sig, Tobs_hr);
    else
        filename = sprintf('SNRgfdotTemplateSpacingvsStepSizefdot_%0.eTobs_%0.f.csv', fdot_sig, Tobs_hr);
    end
end
data = readmatrix(filename);

step_sizes = data(:, 1);
templateSpacings = data(:, 2:end);
gvec = 1:length(templateSpacings(1, :));

% Choose grouping number for plot
g = gvec(15);

%% Plot SNRg Template Spacing vs Step Size 

% Plot 
figure(1)
sk = scatter(step_sizes, templateSpacings(:, g), '*k');
sk.LineWidth = 3;
title(['SNRg_{', num2str(g), '} Template Spacing vs Step Size']);
if (strcmp(TS, 'f'))
    xlabel('f Step Size (Hz)');
    ylabel('f Template Spacing (Hz)')
else
    xlabel('fdot Step Size (Hz/s)');
    ylabel('fdot Template Spacing (Hz/s)');
end
grid on;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 15;