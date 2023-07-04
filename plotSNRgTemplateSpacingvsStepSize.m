%% Plotting SNRg Template Spacing vs fdot Step Size used to calculate template spacing

clear;
clc;
close all;

%% Read Data
% Choose frequency derivative (Hz/s), observation time (hr)
fdot_sig = -5.e-9;
Tobs_hr = 39.;

% Read in data and assign to variables
filename = sprintf('SNRgTemplateSpacingvsStepSizefdot_%0.eTobs_%0.fP.csv', fdot_sig, Tobs_hr);
data = readmatrix(filename);

fdot_sigStepSizes = data(:, 1);
templateSpacings = data(:, 2:end);
ivec = 1:length(templateSpacings(1, :));

% Choose grouping number for plot
i = ivec(end); %3.125e-11 to end

%% Plot SNRg Template Spacing vs Step Size 

% Plot 
figure(1)
sk = scatter(fdot_sigStepSizes, templateSpacings(:, i), '*k');
sk.LineWidth = 3;
title(['SNRg_{', num2str(i), '} Template Spacing vs fdot Step Size']);
xlabel('fdot Step Size (Hz/s)');
ylabel('Template Spacing (Hz/s)');
grid on;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 15;