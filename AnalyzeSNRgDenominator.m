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
Data = Data(1:50, :);
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
noiseMeanDS = noise_mean_DS(noiseMeanPower, 1, Nseg, Nseg);
noiseStdDS = noise_std_DS(noiseMeanPower, 1, Nseg, Nseg);
ivec = 1:Nseg;
iVals = 1:Nseg;
perfiVals = iVals(Nseg./iVals == floor(Nseg./iVals));
preddenoms = noiseStdDS;
perfpreddenoms = sqrt(perfiVals).*noiseMeanPower*sqrt(Nseg);

% Plot
figure;
hold on;
s1 = scatter(perfiVals, perfpreddenoms, 100, 'p');
p1 = plot(ivec, preddenoms, 'LineWidth', 3);
s1.CData = 0.8*[204, 204, 255]/255;
p1.Color = s1.CData;
e1 = errorbar(iVals, denoms, sigmadenoms, sigmadenoms);
e1.LineStyle = 'none';
e1.LineWidth = 2;
title('SNR$_g$ Denominator vs $g$', 'Interpreter', 'LaTex');
xlabel('$g$', 'Interpreter', 'LaTex');
ylabel('SNR$_g$ Denominator', 'Interpreter', 'LaTex');
legend('Predicted $g^*$', 'Predicted', 'Measured', 'Location', 'Northwest', 'Interpreter', 'LaTex');
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 22;
grid on;
hold off;


%% Functions

% Calculate root mean square h_0 for pure signal spectrogram
function [h_0] = h0rms(S, fIndex, Nseg, nBinSide)
    h_0sqr = 0;
    for tIndex = 1:Nseg
        h_0sqr = h_0sqr + sum(abs(S((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2);
    end
    h_0 = sqrt(h_0sqr/Nseg);
end

% Calculate excess segments
function [ib] = ibar(i, Nseg)
    ib = Nseg - floor(Nseg./i).*i;
end

% Calculate mean noise detection statistic from weights and signal spectrogram
function [D_n] = noise_mean_DS(A, M, N, Nseg)
    i = M:N;
    D_n = (floor(Nseg./i) + ibar(i, Nseg)).*A;
end

% Calculate standard deviation of noise detection statistic from weights and signal spectrogram
function [sigma_n] = noise_std_DS(A, M, N, Nseg)
    i = M:N;
    var_n = (floor(Nseg./i).*i.^2 + ibar(i, Nseg).^2).*A^2;
    sigma_n = sqrt(var_n);
end

% Calculate standard deviation of noise detection statistic from weights and signal spectrogram
function [D_s] = pred_signal_DS(h_0, A, M, N, Nseg)
    i = M:N;
    D_s = (floor(Nseg./i).*i.^2 + ibar(i, Nseg).^2).*h_0^2 + noise_mean_DS(A, M, N, Nseg);
end
