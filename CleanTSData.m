%% Script to remove i = 1 template spacings from csv files
% DO NOT RUN THIS SCRIPT UNLESS YOU NEED TO PERMANENTLY DELETE TEMPLATE SPACINGS FROM YOUR CSV FILES
clear;
clc;
close all;

ParsevalSNR = 1;
TS = 'f'; % Choose 'f' or 'fdot'


%% Clean data

% Set fdot(Hz/s)
fdot_sig = -5.e-5;

% Read in data
if (ParsevalSNR == 0)
    filename = sprintf('SNRg%sTemplateSpacingfdot_%0.e.csv', TS, fdot_sig);
else
    filename = sprintf('SNRg%sTemplateSpacingfdot_%0.eP.csv', TS, fdot_sig);
end
Data = readmatrix(filename);

% Find i = 1 rows
iIs1 = (Data(:, 2) == 1);

% Remove rows
Data = Data(~iIs1, :);

% Rewrite file
writematrix(Data, filename);
fprintf('%d Data Points Deleted\n', sum(iIs1));