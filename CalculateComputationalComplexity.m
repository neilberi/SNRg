% Calculates computational complexity of semi-coherent search over fdot and
% Tobs with segment grouping

clear;
clc;
close all;

%% Calculate template count

% Choose grouping number i
i = 10;

% Choose min and max observation periods (in hours)
T1 = 4;
T2 = 32;

% Choose min and max frequency derivatives (Hz/s)
% NOTE: fdot1 < fdot2, but |fdot1| > |fdot2| as fdot < 0
fdot1 = -5.e-5;
fdot2 = -5.e-9;

% Load parameters and calculate template count from analytically integated
% formula
load TemplateSpacingModelParams.mat


