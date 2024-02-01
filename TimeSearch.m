%% Script to time search and compare to calculation from model

%clear;
%clc;
%close all;

ParsevalSNR = 0;

%% Initialize search

% Define search domain
f1 = 10;
f2 = 10.0005;
fdot1 = -1.e-8;
fdot2 = -5.e-9;

% Set observation period
Tobs_hr = 10;
Tobs = 3600*Tobs_hr;

% Set maximum group length
gmax = 5;

% Load search parameters and define model functions
model_func = @(params, Tobs_hr, g, fdot) 10.^( params(1).*log10(Tobs_hr) + ...
                                               params(2).*log10(g) + ... 
                                               params(3).*log10(abs(fdot)) + ... 
                                               params(4) );

model_funcng = @(params, fdot, Tobs_hr) 10^( params(1).*log10(fdot) + ...
                                             params(2).*log10(Tobs_hr) + ... 
                                             params(3) );

params = readmatrix('SNRgModelParams.csv');
params_TSf = params(1, 1:4);
params_TSfdot = params(2, 1:4);
params_SNR = params(3, 1:4);
params_comptime = params(4, 1:3);

% Define grid of templates
Tfdot = fdot2;
while (Tfdot(end) > fdot1)
    Tfdot = [Tfdot, Tfdot(end) - model_func(params_TSfdot, Tobs_hr, min([gmax, floor(Tobs*sqrt(2*abs(Tfdot(end))))]), Tfdot(end))]; %#ok<*AGROW>
end
Tf = cell(1, length(Tfdot));
SNRg = Tf;
for k = 1:length(Tf)
    Tf{k} = f1;
    while (Tf{k}(end) < f2)
        Tf{k} = [Tf{k}, Tf{k}(end) + model_func(params_TSf, Tobs_hr, min([gmax, floor(Tobs*sqrt(2*abs(Tfdot(k))))]), Tfdot(k))]; %#ok<*AGROW>
    end
    SNRg{k} = zeros(size(Tf{k}));
end

% Count templates
N_T = 0;
for k = 1:length(Tfdot)
    N_T = N_T + length(Tf{k});
end
fprintf('For f on [%1.1e, %1.1e] Hz and fdot on [%1.1e, %1.1e] Hz/s,\n', f1, f2, fdot1, fdot2);
fprintf('%e templates are needed.\n\n', N_T);

 % Set information on noise samples
nNoiseSample = 10;
noiseOffset = 100;
nBinSide = 4;

% Sampling frequency
fsamp = 32.;

% Noise ASD ( / sqrt(Hz) )
hnoise = 12;
noiseamp = hnoise;

% Signal amplitude (no antenna pattern) 
hamp = 1;

% Signal initial frequency (Hz)
f_sig = 10+0.005/2.;

% Signal frequency derivative (Hz/s)
fdot_sig = -2.4e-8;
fprintf('Signal parameters: f = %e Hz, fdot = %e Hz/s\n\n', f_sig, fdot_sig);

% Define time series to hold raw data stream of signal plus noise
deltat = 1./fsamp;
t = 0:deltat:Tobs;
Nsample = length(t);

 % Generate noise
noise = noiseamp*random('norm',0.,1.,1.,Nsample);

% Generate signal in time domain
signal = hamp*sin(2*pi*(f_sig*t + 0.5*fdot_sig*t.^2));
data = signal + noise;

%% Iterate through each template and calculate SNRg
timeitT_comp = 0;
%profile on -timer 'real';
for i = 1:length(Tfdot)
    fdot = Tfdot(i);

    % Calculate coherence time 
    Tcoh = 1/sqrt(2*abs(Tfdot(i)));
    Tcoh = floor(Tcoh);
    Tcoh_hr = Tcoh/3600;

    % Calculate Nseg and round Tobs
    Nseg = floor(Tobs/Tcoh);
    Tobs_r = Nseg*Tcoh;
    Tobs_hr_r = Tcoh*Nseg;
    Nsample_r = Tobs_r*fsamp;

    % Segment and transform data
    Nseg = floor(Tobs_r/Tcoh);
    Nsample_coh = floor(Nsample_r/ Nseg);
    rawfft = zeros(Nseg, Nsample_coh);
    for seg = 1:Nseg
       indlo = (seg-1)*Nsample_coh + 1;
       indhi = indlo + Nsample_coh - 1;
       segment = data(indlo:indhi);
       rawfft(seg, :) = fft(segment,Nsample_coh);
    end

    % Set noise power mean from Parseval's Theorem
    noiseMeanPower = Tcoh*fsamp*noiseamp^2;

    % Calculate predicted nosie DS
    if (ParsevalSNR == 1)
        noiseMeanDS = noise_mean_DS(noiseMeanPower, 1, Nseg, Nseg);
        noiseStdDS = noise_std_DS(noiseMeanPower, 1, Nseg, Nseg);
    else
        noiseMeanDS = [];
        noiseStdDS = [];
    end

    for j = 1:length(Tf{i})
        f = Tf{i}(j);

        % Search band width (Hz)
        bandscale = 130; % Increase if fIndex is too high/low
        flo = min(f,f+fdot*Tobs_r);
        fhi = max(f,f+fdot*Tobs_r);
        freqlo_approx = f+fdot*Tobs_r - (noiseOffset + nNoiseSample*(2*nBinSide+1))/Tcoh - bandscale/Tcoh;
        freqhi_approx = f + (noiseOffset + nNoiseSample*(2*nBinSide+1))/Tcoh + bandscale/Tcoh;
        freqhi = ceil(freqhi_approx * Tcoh)/Tcoh;
        freqlo = floor(freqlo_approx * Tcoh)/Tcoh;
        bandwidth = freqhi - freqlo;

        % Generate spectra for each coherence time and extract band of interest to make spectrogram
        indbandlo = floor(freqlo*Tcoh)+1;
        indbandhi = floor(freqhi*Tcoh)+1;
        nbandbin = indbandhi - indbandlo + 1;
        spectrogram = zeros(nbandbin,Nseg);
        rawSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
        for seg = 1:Nseg
           spectrogram(1:nbandbin,seg) = abs(rawfft(seg, indbandlo:indbandhi));
           rawSpectrogram(1:nbandbin,seg) = rawfft(seg, indbandlo:indbandhi);
        end

        % Create vector of indices for frequency trajectory in the middle of each time segment
        freqTrajIndex = round(((1:Nseg) - 0.5).*Nsample_coh) + 1;
        tMid = t(freqTrajIndex);
        
        SNRg{i}(j) = GenerateTemplate(fdot, min([gmax, floor( Tobs_r*sqrt(2*abs(fdot)) )]), rawSpectrogram, f, freqlo, freqhi, hamp, Tcoh, t, tMid, Nseg, nbandbin, indbandlo, indbandhi, Nsample_coh, nNoiseSample, nBinSide, noiseOffset, noiseMeanDS, noiseStdDS, ParsevalSNR);
        %func = @ () GenerateTemplate(fdot, min([gmax, floor( Tobs_r*sqrt(2*abs(fdot)) )]), rawSpectrogram, f, freqlo, freqhi, hamp, Tcoh, t, tMid, Nseg, nbandbin, indbandlo, indbandhi, Nsample_coh, nNoiseSample, nBinSide, noiseOffset, noiseMeanDS, noiseStdDS, ParsevalSNR);
        %timeitT_comp = timeitT_comp + timeit(func);
    end
end

% p = profile('info');
% T_comp = 0;
% for x = 1:length(p.FunctionTable)
%    T_comp = T_comp + p.FunctionTable(x).TotalTime;
% end
% fprintf('Total time (profiler) = %e s\n', T_comp);
% fprintf('Total time (timeit + ignore initialization) = %e s\n', timeitT_comp);
% 
% 
% %% Plot
% figure;
% hold on;
% for i = 1:length(SNRg)
%     scatter3(Tfdot(i).*ones(size(Tf{i})), Tf{i}, SNRg{i});
% end
% 
% %%
% maxes = zeros(size(Tfdot));
% for i = 1:length(maxes)
%     maxes(i) = max(SNRg{i});
% end
% [x, fdotIndexMax] = max(maxes);
% [maxSNR, fIndexMax] = max(SNRg{fdotIndexMax});



%% Functions

% Calculate excess segments
function [gb] = gbar(g, Nseg)
    gb = Nseg - floor(Nseg./g).*g;
end

% Calculate mean noise detection statistic from weights and signal spectrogram
function [D_n] = noise_mean_DS(A, M, N, Nseg)
    g = M:N;
    D_n = (floor(Nseg./g) + gbar(g, Nseg)).*A;
end

% Calculate standard deviation of noise detection statistic from weights and signal spectrogram
function [sigma_n] = noise_std_DS(A, M, N, Nseg)
    g = M:N;
    var_n = (floor(Nseg./g).*g.^2 + gbar(g, Nseg).^2).*A^2;
    sigma_n = sqrt(var_n);
end

% Calculate standard deviation of noise detection statistic from weights and signal spectrogram
function [D_s] = pred_signal_DS(h_0, A, M, N, Nseg)
    g = M:N;
    D_s = (floor(Nseg./g).*g.^2 + gbar(g, Nseg).^2).*h_0^2 + noise_mean_DS(A, M, N, Nseg);
end
                            