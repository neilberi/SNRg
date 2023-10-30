%% Script to Plot and Analyze the Variance of the Denominator of SNRg_i

clear;
clc;
close all;

ParsevalSNR = 0;

%% Read in Data

% Set signal amplitude and noise ASD in sqrt(1/Hz)
hamp = 1.e-3;
hnoise = 100;
noiseamp = hnoise;

% Set fsamp (Hz), fdot (Hz/s), and Tobs (hr)
fsamp = 32;
fdot_sig = -5.e-5;
Tobs_hr = 16;
Tobs = Tobs_hr*3600;

% Read in data from csv
if (ParsevalSNR == 1)
    filename4 = sprintf('SNRgDenominatorfdot_%0.es-2Tobs_%0.fhrnoiseamp_%.2esqrt(s)sigamp_%.2eP.csv', fdot_sig, Tobs_hr, hnoise, hamp);
else
    filename4 = sprintf('SNRgDenominatorfdot_%0.es-2Tobs_%0.fhrnoiseamp_%.2esqrt(s)sigamp_%.2e.csv', fdot_sig, Tobs_hr, hnoise, hamp);
end
Data = readmatrix(filename4);
Ntrials = height(Data);

% Calculate and round Tcoh (s)
if (abs(fdot_sig) > 0)
Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
else 
    Tcoh_hr = 1./3600.*sqrt(0.5/abs(-5.e-6));
end
Tcoh = Tcoh_hr * 3600.;
Tcoh = floor(Tcoh);
Tcoh_hr = Tcoh/3600.;
Nseg = floor(Tobs/Tcoh);
Tobs = Nseg*Tcoh;
Tobs_hr = Tobs / 3600.;

% Calculate means and standard deviations for each i
denoms = mean(Data);
sigmadenoms = std(Data)/sqrt(Ntrials);

% Calculate predicted values
noiseMeanPower = Tcoh*fsamp*noiseamp^2;
noiseStdDS = noiseMeanPower*sqrt(Nseg);
finei = linspace(1, Nseg, 1000);
iVals = 1:Nseg;
perfiVals = iVals(Nseg./iVals == floor(Nseg./iVals));
preddenoms = sqrt(finei).*noiseStdDS;
perfpreddenoms = sqrt(perfiVals).*noiseStdDS;

% Plot
figure;
hold on;
s1 = scatter(perfiVals, perfpreddenoms, 80, 'p');
p1 = plot(finei, preddenoms);
s1.CData = 0.8*[204, 204, 255]/255;
p1.Color = s1.CData;
e1 = errorbar(iVals, denoms, sigmadenoms, sigmadenoms);
e1.LineStyle = 'none';
e1.LineWidth = 2;
title('SNRg_i Denominators vs i');
xlabel('i');
ylabel('SNR Denominator');
legend('Predicted (i divides N_{seg})', 'Predicted', 'Measured', 'Location', 'Northwest');
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 16;
grid on;
hold off;


