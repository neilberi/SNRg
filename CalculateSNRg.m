% SNRg Calculations

clear;
clc;
close all;

% Choose settings (0 = off, 1 = on)
ParsevalSNR = 1; 
RealPartDS = 0;
Animation = 0;
NequalNseg = 0;
PlotSNRvsdfdot = 1;
PlotMaxSNRgvsTobs = 0;
PlotLastTimeSegment = 0;
Plotdfdotvsi = 0;
OutputFile1 = 0;
OutputFile2 = 0;
OutputFile3 = 0;

%% Generate Signal Spectrogram

% SNR1 and {SNRg_i | M<=i<=N} will be calculated
M = 47;
N = 47;
Ntrial = 10;
if (OutputFile1 == 0)
    Ntrial = 1;
end

% Set information on noise samples
nNoiseSample = 10;
noiseOffset = 50;
nBinSide = 4;

% Sampling frequency
fsamp = 1024.;

% Noise ASD ( / sqrt(Hz) )
hnoise = 12;
noiseamp = hnoise;

% Signal initial frequency (Hz)
f_sig = 100.;

% Signal frequency derivative (Hz/s)
fdot_sig = -5.e-7;

% Length of observation (hr)
Tobs_hr = 25.;
Tobs = Tobs_hr * 3600.;
if (OutputFile1 == 1)
    if (ParsevalSNR == 1)
        filename = sprintf('SNRgfdot_%0.es-2Tobs_%0.fhrP.csv', fdot_sig, Tobs_hr);
    else
        filename = sprintf('SNRgfdot_%0.es-2Tobs_%0.fhr.csv', fdot_sig, Tobs_hr);
    end
end

% Coherence time (hr) - choose so that signal drifts 0.5 bin per coherence time
Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
Tcoh = Tcoh_hr * 3600.;

% Adjust coherence time to be an integer number of seconds
Tcoh = floor(Tcoh);
Tcoh_hr = Tcoh/3600.;

% DO NOT TOUCH/CARE
searchScale = 0;

f_sigStepSize = (0.1)/Tcoh*(Tcoh < Tobs) + 0.3/sqrt(0.5/(1.e-6))*(Tcoh >= Tobs);
fvec_sig = (f_sig-searchScale*f_sigStepSize) : f_sigStepSize : (f_sig+searchScale*f_sigStepSize);

fdot_sigStepSize = (7.8e-6)*abs(fdot_sig)*(fdot_sig ~= 0) + 1.e-7*(fdot_sig == 0);
fdotvec_sig = (fdot_sig-searchScale*fdot_sigStepSize) : fdot_sigStepSize : (fdot_sig+searchScale*fdot_sigStepSize);

[F, Fdot] = meshgrid(fvec_sig, fdotvec_sig);

% Adjust observation time to be exactly an integer number of coherence times to avoid binning headaches
Nseg = floor(Tobs/Tcoh);
if (NequalNseg == 1)
    N = Nseg;
    %N = floor(Nseg/10)*10;
end
Tobs = Nseg*Tcoh;
Tobs_hr = Tobs / 3600.;
fprintf('Observation time = %d sec (%f hr)\n',Tobs,Tobs_hr);
fprintf('Coherence time per segment = %d sec (%f hr)\n',Tcoh,Tcoh_hr);
assert(M<=N, 'Error: M>N');

% Set noise power mean and standard deviation from Parseval's Theorem
if (ParsevalSNR == 1)
    noiseMean = 1;
    noiseSTD = 1;
end

% Define time series to hold raw data stream of signal plus noise
deltat = 1./fsamp;
t = 0:deltat:Tobs;
Nsample = length(t);

% Search band width (Hz)
bandscale = 115; % Increase if fIndex is too high/low
f_siglo = min(f_sig,f_sig+fdot_sig*Tobs);
f_sighi = max(f_sig,f_sig+fdot_sig*Tobs);
freqlo_approx = min(fvec_sig)+min(fdotvec_sig)*Tobs - (noiseOffset + nNoiseSample + nBinSide)/Tcoh - bandscale/Tcoh;
freqhi_approx = max(fvec_sig)+max(fdotvec_sig)*Tobs*(max(fdotvec_sig) > 0) + (noiseOffset + nNoiseSample + nBinSide)/Tcoh + bandscale/Tcoh;
freqhi = ceil(freqhi_approx * Tcoh)/Tcoh;
freqlo = floor(freqlo_approx * Tcoh)/Tcoh;
bandwidth = freqhi - freqlo;
fprintf('Low/high frequencies of band shown: %f-%f Hz\n',freqlo,freqhi);

% Number of bins signal drifts per coherence time
Nbin_drift = abs(fdot_sig) * Tcoh * Tcoh;
Nbin_drift_total = abs(fdot_sig) * Tobs * Tcoh;
fprintf('Number of bins signal drifts per coherence time = %8.2f\n',Nbin_drift);
fprintf('Total number of bins signal drifts = %8.2f\n',Nbin_drift_total);

for dummyIterator = 1:Ntrial
    % Generate noise
    noise = noiseamp*random('norm',0.,1.,1.,Nsample);

    % Signal amplitude (no antenna pattern) 
    hamp = 1;

    % Generate signal in time domain
    signal = hamp*sin(2*pi*(f_sig*t + 0.5*fdot_sig*t.^2));
    data = signal + noise;

    % Generate spectra for each coherence time and extract band of interest to make spectrogram
    indbandlo = floor(freqlo*Tcoh)+1;
    indbandhi = floor(freqhi*Tcoh)+1;
    nbandbin = indbandhi - indbandlo + 1;
    Nseg = floor(Tobs/Tcoh);
    Nsample_coh = floor(Nsample / Nseg);
    spectrogram = zeros(nbandbin,Nseg);
    rawSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
    fftNormalization = sqrt(1/Tcoh/fsamp)/noiseamp;
    for seg = 1:Nseg
       fprintf('Generating segment %d\n',seg);
       indlo = (seg-1)*Nsample_coh + 1;
       indhi = indlo + Nsample_coh - 1;
       segment = data(indlo:indhi);
       segment_noNoise = signal(indlo:indhi);
       rawfft = fft(segment,Nsample_coh);
       spectrogram(1:nbandbin,seg) = abs(rawfft(indbandlo:indbandhi))*fftNormalization;
       rawSpectrogram(1:nbandbin,seg) = rawfft(indbandlo:indbandhi)*fftNormalization;
       [dummyVariable, fIndex_sig] = max(rawfft(indbandlo:indbandhi)*fftNormalization);
       if (seg == round(Nseg/5) || seg == round(Nseg/5*2) || seg == round(Nseg/5*3) || seg == round(Nseg/5*4) || seg == Nseg)
           fprintf('Signal fIndex of segment %d = %d; magnitude = %d; phase = %d\n', seg, fIndex_sig, abs(dummyVariable), angle(dummyVariable));
       end
    end

    % Create modified spectrogram for plot
    plotSpectrogram = spectrogram;
    for k = 1:Nseg
        segment = plotSpectrogram(:, k);
        segment(segment < 0.1*max(segment)) = 0;
        plotSpectrogram(:, k) = segment;
    end

    % Plot spectrogram bar3(spectrogram)
    segarray = 1:Nseg;
    seghour = (segarray-1)*Tcoh/3600.;
    indarray = [indbandlo:indbandhi];
    freqplot = (indarray-indbandlo)*1.0/Tcoh + freqlo;

    figure(1);
    fprintf('Creating spectrogram plot with %d time bins (columns) and %d frequency bins (rows)...\n',Nseg,nbandbin);
    s1 = surf(seghour,freqplot,plotSpectrogram);
    title('Raw spectrogram');
    xlabel('Time (hours)');
    ylabel('Frequency (Hz)');
    view(0,90);
    ax = gca;
    ax.FontSize = 18;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);

    %% Calculate Detection Statistics

    % Create vector of column midpoints
    tMid = zeros(1, Nseg);
    for tIndex = 1:Nseg
    tMid(tIndex) = (tIndex - 0.5)*Tcoh;
    end

    tMid_hr = tMid./3600;

    % Create vector of indices for frequency trajectory
    freqTrajIndex = zeros(1, Nseg);
    for tIndex = 1:Nseg
        freqTrajIndex(tIndex) = find(t == tMid(tIndex));
    end

    % Calculate expected phase/frequency trajectories and signal
    freqTraj = f_sig + fdot_sig.*t;
    phaseTraj = 2*pi.*(f_sig.*t + 0.5*fdot_sig.*t.^2);
    searchSignal = hamp.*sin(phaseTraj);

    % Create template spectrogram
    searchSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
    for seg = 1:Nseg
        indlo = (seg-1)*Nsample_coh + 1;
        indhi = indlo + Nsample_coh - 1;
        searchSegment = searchSignal(indlo:indhi);
        searchRawfft = fft(searchSegment, Nsample_coh);
        searchSpectrogram(1:nbandbin,seg) = searchRawfft(indbandlo:indbandhi)*fftNormalization;
    end

    % Find predicted frequency at each column midpoint
    fMid = zeros(1, Nseg);
    for tIndex = 1:Nseg
        fMid(tIndex) = freqTraj(freqTrajIndex(tIndex));
    end

    % Calculate signal frequency indices
    fIndex = zeros(1, Nseg);
    for tIndex = 1:Nseg
        fIndex(tIndex) = find(abs(freqplot - fMid(tIndex)) == min(abs(freqplot - fMid(tIndex))), 1, 'first');
    end

    % Calculate background noise frequency indices
    fIndex_noise = zeros(nNoiseSample, Nseg);
    for k = 1:nNoiseSample
        fIndex_noise(k, :) = fIndex + (2*nBinSide+2)*k + noiseOffset;
    end

    % Calculate weights array
    weight = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
    for tIndex = 1:Nseg
        %weight(:, tIndex) = hamp.*searchSpectrogram(:, tIndex)/sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2);
        weight(:, tIndex) = searchSpectrogram(:, tIndex).*abs(searchSpectrogram(:, tIndex)).^2/sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2);
    end

    % Calculate detection statistics
    signalDS1 = 0;
    signalDSg_i = zeros(N-M+1, 1);
    grouper_sig = (1+sqrt(-1)).*zeros(N-M+1, 1);
    noiseDS1 = zeros(1, nNoiseSample);
    noiseDSg_i = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);
    grouper_noise = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);

    signalDSContribution1 = [];
    signalDSContribution2 = [];
    signalDSContribution3 = [];
    signalDSContribution4 = [];

    noiseDSContribution1 = [];
    noiseDSContribution2 = [];
    noiseDSContribution3 = [];
    noiseDSContribution4 = [];

    for tIndex = 1:Nseg
        % Signal detection statistics
        signalDS1 = signalDS1 + (abs(spectrogram(fIndex(tIndex), tIndex)))^2;
        grouper_sig = grouper_sig + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex));
        for i = M:N
            if (RealPartDS == 1)
                contribution = real(grouper_sig(i-M+1)); 
            else
                contribution = abs(grouper_sig(i-M+1))^2;
            end
            if (mod(tIndex, i) == 0)
                signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + contribution;
                if (i == 1)
                    signalDSContribution1 = [signalDSContribution1, contribution];
                elseif (i == 2)
                    signalDSContribution2 = [signalDSContribution2, contribution];
                elseif (i == 3)
                    signalDSContribution3 = [signalDSContribution3, contribution];
                elseif (i == 4)
                    signalDSContribution4 = [signalDSContribution4, contribution];
                end
                grouper_sig(i-M+1) = 0;
            elseif (tIndex == Nseg)
                signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + contribution;
                grouper_sig(i-M+1) = 0;
            end
        end

        % Noise detection statistics
        for k = 1:nNoiseSample
            noiseDS1(k) = noiseDS1(k) + (abs(spectrogram(fIndex_noise(k, tIndex), tIndex)))^2;
            grouper_noise(:, k) = grouper_noise(:, k) + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex_noise(k, tIndex)-nBinSide):(fIndex_noise(k, tIndex)+nBinSide), tIndex));
        end
        for i = M:N
            if (RealPartDS == 1)
                contribution = real(grouper_noise(i-M+1, :)); 
            else
                contribution = abs(grouper_noise(i-M+1, :)).^2;
            end
            if (mod(tIndex, i) == 0)
                noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + contribution;
                if (i == 1)
                    noiseDSContribution1 = [noiseDSContribution1, [abs(grouper_noise(i-M+1, 1))^2; abs(grouper_noise(i-M+1, floor(nNoiseSample/2)))^2; abs(grouper_noise(i-M+1, end))^2]];
                elseif (i == 2)
                    noiseDSContribution2 = [noiseDSContribution2, [abs(grouper_noise(i-M+1, 1))^2; abs(grouper_noise(i-M+1, floor(nNoiseSample/2)))^2; abs(grouper_noise(i-M+1, end))^2]];
                elseif (i == 3)
                    noiseDSContribution3 = [noiseDSContribution3, [abs(grouper_noise(i-M+1, 1))^2; abs(grouper_noise(i-M+1, floor(nNoiseSample/2)))^2; abs(grouper_noise(i-M+1, end))^2]];
                elseif (i == 4)
                    noiseDSContribution4 = [noiseDSContribution4, [abs(grouper_noise(i-M+1, 1))^2; abs(grouper_noise(i-M+1, floor(nNoiseSample/2)))^2; abs(grouper_noise(i-M+1, end))^2]];
                end
                grouper_noise(i-M+1, :) = 0;
            elseif (tIndex == Nseg)
                noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + contribution;
            end
        end
    end

    %% Calculate SNRs and Plot
    % Calculate SNRs
    SNR1 = abs(signalDS1 - mean(noiseDS1))/std(noiseDS1);
    fprintf('SNR1 = %e\n', SNR1);
    SNRg_i = zeros(N-M+1, 1);
    for i = M:N
        if (ParsevalSNR == 0)
            SNRg_i(i-M+1) = abs(signalDSg_i(i-M+1) - mean(noiseDSg_i(i-M+1, :)))/std(noiseDSg_i(i-M+1, :));
        elseif (ParsevalSNR == 1)
            SNRg_i(i-M+1) = abs(signalDSg_i(i-M+1) - noiseMean)/noiseSTD/sqrt(i);
        end
        fprintf('SNRg_%d = %e\n', i, SNRg_i(i-M+1));
    end

    if (OutputFile1 == 1)
        assert(all(M:N==1:Nseg), 'Error: Output File 1 can only function when M=1 and N=Nseg')
        fout = fopen(filename, 'a');
        for i = M:N
            if (i == N)
                fprintf(fout, '%e\n', SNRg_i(i-M+1));
            else
                fprintf(fout, '%e,', SNRg_i(i-M+1));
            end
        end
    end

end
% Plot data
figure(2)
p1 = scatter(M:N, SNRg_i, 'ok');
hold on;
if (ParsevalSNR == 0)
    p2 = plot(M:N, SNR1.*ones(1, N-M+1), '--r');
end
title('SNRg_i vs i');
xlabel('i');
ylabel('SNRg_i');
if (ParsevalSNR == 0)
    legend('SNRg_i', 'SNR1', 'location', 'northwest');
else
    legend('SNRg_i', 'location', 'northwest');
end
grid on;
p1.LineWidth = 3;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 16;

[imax, SNRg_imax] = max(SNRg_i);
fprintf('Max of %e at i = %d\n', imax, SNRg_imax);
%% Plot DS Contributions

if (NequalNseg == 1 && M == 1 && Nseg >= 4)
    % Create segment vectors for each SNRg_i
    segvecs = {[], [], [], []};
    for tIndex = 1:Nseg
        for i = 1:4
            if (mod(tIndex, i) == 0)
                segvecs{i} = [segvecs{i}, tIndex];
            end
        end
    end

    figure(3)
    subplot(2, 2, 1)
    p311 = plot(segvecs{1}, signalDSContribution1, '-r');
    hold on;
    p321 = plot(segvecs{2}, signalDSContribution2, '-g');
    p331 = plot(segvecs{3}, signalDSContribution3, '-b');
    p341 = plot(segvecs{4}, signalDSContribution4, '-m');
    hold off;
    ylim([0, 17]);
    title('Signal DS Contributions for SNRg_i');
    xlabel('time segment');
    ylabel('contribution');
    legend('i=1', 'i=2', 'i=3', 'i=4');
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 15;

    subplot(2, 2, 2)
    p312 = plot(segvecs{1}, noiseDSContribution1(1, :), '-r');
    hold on;
    p322 = plot(segvecs{2}, noiseDSContribution2(1, :), '-g');
    p332 = plot(segvecs{3}, noiseDSContribution3(1, :), '-b');
    p342 = plot(segvecs{4}, noiseDSContribution4(1, :), '-m');
    hold off;
    title('Noise DS Contributions for SNRg_i (first noise sample)');
    xlabel('time segment');
    ylabel('contribution');
    legend('i=1', 'i=2', 'i=3', 'i=4');
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 15;

    subplot(2, 2, 3)
    p313 = plot(segvecs{1}, noiseDSContribution1(2, :), '-r');
    hold on;
    p323 = plot(segvecs{2}, noiseDSContribution2(2, :), '-g');
    p333 = plot(segvecs{3}, noiseDSContribution3(2, :), '-b');
    p343 = plot(segvecs{4}, noiseDSContribution4(2, :), '-m');
    hold off;
    title('Noise DS Contributions for SNRg_i (middle noise sample)');
    xlabel('time segment');
    ylabel('contribution');
    legend('i=1', 'i=2', 'i=3', 'i=4');
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 15;

    subplot(2, 2, 4)
    p314 = plot(segvecs{1}, noiseDSContribution1(3, :), '-r');
    hold on;
    p324 = plot(segvecs{2}, noiseDSContribution2(3, :), '-g');
    p334 = plot(segvecs{3}, noiseDSContribution3(3, :), '-b');
    p344 = plot(segvecs{4}, noiseDSContribution4(3, :), '-m');
    hold off;
    title('Noise DS Contributions for SNRg_i (last noise sample)');
    xlabel('time segment');
    ylabel('contribution');
    legend('i=1', 'i=2', 'i=3', 'i=4');
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 15;

end
%% Calculate SNR for Various fdot

TS_model = @(a, fdot_log, Tobs_log, i_log) 10^( (a(1) + a(2).*fdot_log).*Tobs_log + ...
                                                (a(3) + a(4).*fdot_log).*i_log + ...
                                                (a(5) + a(6).*fdot_log) );
load TemplateSpacingModelParams.mat

% Create vector of fdot values (for search)
searchScale = 5;
fdot_sigStepSize = TS_model(TemplateSpacingModelParams, log10(abs(fdot_sig)), log10(Tobs_hr), log10(N));
fdotvec_sig = (fdot_sig - searchScale*fdot_sigStepSize):fdot_sigStepSize:(fdot_sig + searchScale*fdot_sigStepSize);
SNR1Array = zeros(1, length(fdotvec_sig));
SNRg_iArray = zeros(N-M+1, length(fdotvec_sig));

for j = 1:length(fdotvec_sig)
    
    fprintf('Calulating SNRg_is for dfdot = %e\n', fdotvec_sig(j)-fdot_sig);
    
    % Calculate expected phase/frequency trajectories and signal
    freqTraj = f_sig + fdotvec_sig(j).*t; 
    phaseTraj = 2*pi.*(f_sig.*t + 0.5*fdotvec_sig(j).*t.^2);
    searchSignal = hamp.*sin(phaseTraj);

    % Create template spectrogram
    searchSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
    for seg = 1:Nseg
        indlo = (seg-1)*Nsample_coh + 1;
        indhi = indlo + Nsample_coh - 1;
        searchSegment = searchSignal(indlo:indhi);
        searchRawfft = fft(searchSegment, Nsample_coh);
        searchSpectrogram(1:nbandbin,seg) = searchRawfft(indbandlo:indbandhi)*fftNormalization;
    end

    % Find predicted frequency at each column midpoint
    fMid = zeros(1, Nseg);
    for tIndex = 1:Nseg
        fMid(tIndex) = freqTraj(freqTrajIndex(tIndex));
    end

    % Calculate signal frequency indices
    fIndex = zeros(1, Nseg);
    for tIndex = 1:Nseg
        fIndex(tIndex) = find(abs(freqplot - fMid(tIndex)) == min(abs(freqplot - fMid(tIndex))), 1, 'first');
    end
    
    % Calculate background noise frequency indices
    fIndex_noise = zeros(nNoiseSample, Nseg);
    for k = 1:nNoiseSample
        fIndex_noise(k, :) = fIndex + 10*k + noiseOffset;
    end

    % Calculate weights array
    weight = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
    for tIndex = 1:Nseg
        weight(:, tIndex) = hamp.*searchSpectrogram(:, tIndex)/sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2);
        weight(:, tIndex) = searchSpectrogram(:, tIndex).*abs(searchSpectrogram(:, tIndex)).^2/sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2);
    end

    % Calculate detection statistics
    signalDS1 = 0;
    signalDSg_i = zeros(N-M+1, 1);
    grouper_sig = (1+sqrt(-1)).*zeros(N-M+1, 1);
    noiseDS1 = zeros(1, nNoiseSample);
    noiseDSg_i = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);
    grouper_noise = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);

    if (j == 15)
        DScontributions15 = zeros(1, Nseg);
    elseif (j == 21)
        DScontributions21 = zeros(1, Nseg);
    elseif (j == 11)
        DScontribution11 = zeros(1, Nseg);
    end

    for tIndex = 1:Nseg
        % Signal detection statistics
        signalDS1 = signalDS1 + (abs(spectrogram(fIndex(tIndex), tIndex)))^2;
        grouper_sig = grouper_sig + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex));
        for i = M:N
            if (RealPartDS == 1)
                contribution = real(grouper_sig(i-M+1)); 
            else
                contribution = abs(grouper_sig(i-M+1))^2;
            end
            if (mod(tIndex, i) == 0)
                signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + contribution;
                if (j == 15 && i == 1)
                    DScontributions15(tIndex) = contribution;
                elseif (j == 21 && i == 1)
                    DScontributions21(tIndex) = contribution;
                elseif (j == 11 && i == 1)
                    DScontribution11(tIndex) = contribution;
                end
                grouper_sig(i-M+1) = 0;
            elseif (tIndex == Nseg)
                signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + contribution;
            end
        end

        % Noise detection statistics
        for k = 1:nNoiseSample
            noiseDS1(k) = noiseDS1(k) + (abs(spectrogram(fIndex_noise(k, tIndex), tIndex)))^2;
            grouper_noise(:, k) = grouper_noise(:, k) + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex_noise(k, tIndex)-nBinSide):(fIndex_noise(k, tIndex)+nBinSide), tIndex));
        end
        for i = M:N
            if (RealPartDS == 1)
                contribution = real(grouper_noise(i-M+1, :)); 
            else
                contribution = abs(grouper_noise(i-M+1, :)).^2;
            end
            if (mod(tIndex, i) == 0)
                noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + contribution;
                grouper_noise(i-M+1, :) = 0;
            elseif (tIndex == Nseg)
                noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + contribution;
            end
        end
    end
    
    SNR1Array(j) = abs(signalDS1 - mean(noiseDS1))/std(noiseDS1);
    for i = M:N
        if (ParsevalSNR == 0)
            SNRg_iArray(i-M+1, j) = abs(signalDSg_i(i-M+1) - mean(noiseDSg_i(i-M+1, :)))/std(noiseDSg_i(i-M+1, :));
        elseif (ParsevalSNR == 1)
            SNRg_iArray(i-M+1, j) = abs(signalDSg_i(i-M+1) - noiseMean)/noiseSTD/sqrt(i);
        end
    end
    
    % Compare search spectrogram and signal spectrogram last time segment
    if (PlotLastTimeSegment == 1)
        range = (fIndex(end)-10):(fIndex(end)+10);
        
        figure(1900+j)
        s1900j = scatter(range, abs(searchSpectrogram(range, end)), '*k');
        hold on;
        p1900j = plot(range, abs(rawSpectrogram(range, end)), '-r');
        x1900j = xline(fIndex(end), '--k');
        x1900j.LineWidth = 1;
        s1900j.LineWidth = 3;
        title(['Last Time Segment Fourier Coefficient Magnitude dfdot_{sig} = ', num2str(fdotvec_sig(j) - fdot_sig), ' Hz/s']);
        xlabel('frequency bin');
        ylabel('magnitude');
        legend('template', 'data', 'fIndex');
        ax = gca;
        ax.FontSize = 14;
        ax.LineWidth = 3;
        grid on;
        hold off;
        
        figure(2900+j)
        s2900j = scatter(range, unwrap(atan2(imag(searchSpectrogram(range, end)), real(searchSpectrogram(range, end)))), '*k');
        hold on;
        p2900j = plot(range, unwrap(atan2(imag(rawSpectrogram(range, end)), real(rawSpectrogram(range, end)))), '-r');
        x2900j = xline(fIndex(end), '--k');
        x2900j.LineWidth = 1;
        s2900j.LineWidth = 3;
        title(['Last Time Segment Fourier Coefficient Phase dfdot_{sig} = ', num2str(fdotvec_sig(j) - fdot_sig), ' Hz/s']);
        xlabel('frequency bin');
        ylabel('phase');
        legend('template', 'data', 'fIndex');
        ax = gca;
        ax.FontSize = 14;
        ax.LineWidth = 3;
        grid on;
        hold off;
    end
end
if (searchScale >= 20)
    figure(3000)
    p39001 = plot(1:Nseg, DScontributions15, '-g');
    hold on;
    p39002 = plot(1:Nseg, DScontributions21, '-m');
    p39003 = plot(1:Nseg, DScontribution11, '-b');
    hold off;
    title('Center Bin and Off Bin DS Contributions(i == 1)')
    xlabel('time segment')
    ylabel('DS Contribution')
    legend('off bin with max SNR', 'center bin', 'off bin with lower SNR')
    ax = gca;
    ax.FontSize = 16;
    ax.LineWidth = 3;
    grid on;
end

%% Plot data for deviations in fdot

% Create titles
titleStrings = cell(1, N-M+1);
for i = M:N
    titleStrings{i-M+1} = sprintf('SNRg_{%d} vs dfdot', i);
end
dfdot_sig = fdotvec_sig-fdot_sig;

% Plot SNRg_i vs i
for k = 1:length(dfdot_sig)
    figure(3+k)
    sk1 = scatter(M:N, SNRg_iArray(:, k), 'ok');
    title(['SNRg_i vs i for dfdot = ', num2str(dfdot_sig(k))]);
    xlabel('i');
    ylabel('SNRg_i');
    grid on;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
end

% Plot SNRg_i vs dfdot
if (PlotSNRvsdfdot == 1)
    for i = M:N
        figure(3+length(dfdot_sig)+i-M+1)
        subplot(1, 2, 1)
        hold on;
        si1 = scatter(dfdot_sig, SNRg_iArray(i-M+1, :)/SNRg_iArray(i-M+1, searchScale+1), 'ok');
        title(titleStrings{i-M+1});
        xlabel('dfdot (Hz/s)');
        ylabel('SNR/SNR_m_a_x');
        legend(['i = ', num2str(i)])
        grid on;
        ax = gca;
        ax.LineWidth = 3;
        ax.FontSize = 16;

        subplot(1, 2, 2)
        si2 = scatter(dfdot_sig, SNRg_iArray(i-M+1, :), 'ok');
        title(titleStrings{i-M+1});
        xlabel('dfdot (Hz/s)');
        ylabel('SNR');
        legend(['i = ', num2str(i)])
        grid on;
        ax = gca;
        ax.LineWidth = 3;
        ax.FontSize = 16;
    end
end

%% Animation
if (Animation == 1)
    figure(4+length(dfdot_sig)+N-M+1)
    sgtitle('SNRg_i vs dfdot for Various i', 'FontSize', 20);
    for repeat = 1:3
        for i = M:N
            subplot(1, 2, 1)
            scatter(dfdot_sig, SNRg_iArray(i-M+1, :)/max(SNRg_iArray(i-M+1, :)), 'ok');
            axis([1.01*min(dfdot_sig), 1.01*max(dfdot_sig), 0, 1]);
            xlabel('dfdot (Hz/s)');
            ylabel('SNR/SNR_m_a_x');
            legend(['i = ', num2str(i)]);
            title('Normalized');
            ax = gca;
            ax.LineWidth = 3;
            ax.FontSize = 16;
            grid on;

            subplot(1, 2, 2)
            scatter(dfdot_sig, SNRg_iArray(i-M+1, :), 'ok');
            axis([1.01*min(dfdot_sig), 1.01*max(dfdot_sig), 0, max(SNRg_iArray(:))]);
            xlabel('dfdot (Hz/s)');
            ylabel('SNR');
            legend(['i = ', num2str(i)]);
            title('Not Normalized');
            ax = gca;
            ax.LineWidth = 3;
            ax.FontSize = 16;
            grid on;
            pause(0.5)
        end
    end
end
%% Fit SNRg_i vs dfdot to curve and find dfdot for 20% offset
dfdotp = linspace(min(dfdot_sig), max(dfdot_sig), 1.e4*length(dfdot_sig));
SNRsplines = zeros(N-M+1, length(dfdotp));
templateSpacings = zeros(N-M+1, 1);

for i = M:N
    SNRsplines(i-M+1, :) = spline(dfdot_sig, SNRg_iArray(i-M+1, :), dfdotp);
    templateSpacings(i-M+1) = mode(abs(dfdotp(abs(SNRsplines(i-M+1, :) - 0.8*SNRg_iArray(i-M+1, searchScale+1)) == min(abs(SNRsplines(i-M+1, :) - 0.8*SNRg_iArray(i-M+1, searchScale+1))))));
    
    if (PlotSNRvsdfdot == 1)
        figure(3+length(dfdot_sig)+i-M+1)
        subplot(1, 2, 1)
        hold on;
        pi1 = plot(dfdotp, SNRsplines(i-M+1, :)/SNRg_iArray(i-M+1, searchScale+1));
        legend(['i = ', num2str(i)], 'spline', 'Location', 'south')

        subplot(1, 2, 2)
        hold on;
        pi2 = plot(dfdotp, SNRsplines(i-M+1, :));
        legend(['i = ', num2str(i)], 'spline', 'Location', 'south')
    end
end


for i = M:N
    fprintf('For i = %d, stepping %e in fdot causes a 20%% drop in SNRg_i\n', i, templateSpacings(i-M+1));
end

figure(1000000)
s1000000 = scatter(log10(M:N), log10(templateSpacings));
s1000000.MarkerFaceColor = 'k';
title('Template Spacing vs i');
xlabel('log(i)');
ylabel('log(dfdot)');
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 3;
grid on;


if (OutputFile2 == 1)
    if (ParsevalSNR == 1)
        filename2 = sprintf('SNRgTemplateSpacingfdot_%0.eP.csv', fdot_sig);
    else
        filename2 = sprintf('SNRgTemplateSpacingfdot_%0.e.csv', fdot_sig);
    end
    fout2 = fopen(filename2, 'a');
    for i = M:N
        fprintf(fout2, '%f, %d, %e\n', Tobs_hr, i, templateSpacings(i-M+1));
    end
end

if (OutputFile3 == 1)
    assert(all(M:N==1:Nseg), 'Error: Output File 3 can only function when M=1 and N=Nseg')
    if (ParsevalSNR == 1)
        filename3 = sprintf('SNRgTemplateSpacingvsStepSizefdot_%0.eTobs_%0.fP.csv', fdot_sig, Tobs_hr);
    else
        filename3 = sprintf('SNRgTemplateSpacingvsStepSizefdot_%0.eTobs_%0.f.csv', fdot_sig, Tobs_hr);
    end
    fout3 = fopen(filename3, 'a');
    fprintf(fout3, '%e,', fdot_sigStepSize);
    for i = M:N
        if (i == N)
            fprintf(fout3, '%e\n', templateSpacings(i-M+1));
        else
            fprintf(fout3, '%e,', templateSpacings(i-M+1));
        end
        
    end
end

%% Plot log-log graph
if (Plotdfdotvsi == 1)
    templateFitCoeffs = polyfit(log10(M:N), log10(templateSpacings), 1);

    figure(5+length(dfdot_sig)+N)
    p5fN = plot(log10(M:N), log10(abs(fitCoeffs*sqrt(log(5/4)))), '-k');
    hold on;
    p5fN2 = plot(log10(M:N), -2*(log10(M:N)-log10(N)) + log10(abs(fitCoeffs(N)*sqrt(log(5/4)))), '-r');
    p5fN1 = plot(log10(M:N), -1*(log10(M:N)-log10(N)) + log10(abs(fitCoeffs(N)*sqrt(log(5/4)))), '-g');
    p5fNh = plot(log10(M:N), -0.5*(log10(M:N)-log10(N)) + log10(abs(fitCoeffs(N)*sqrt(log(5/4)))), '-b');
    p5fNf = plot(log10(M:N), -0.25*(log10(M:N)-log10(N)) + log10(abs(fitCoeffs(N)*sqrt(log(5/4)))), '-m');
    p5fNfit = plot(log10(M:N), templateFitCoeffs(1).*log10(M:N) + templateFitCoeffs(2));
    xlim([log10(M), log10(N)]);
    title('dfdot for 20% Offset in SNRg_i vs i (log-log)');
    xlabel('log(i)');
    ylabel('log(dfdot)');
    legend('dfdot values', 'line through last point with slope = -2', 'line through last point with slope = -1', 'line through last point with slope = -1/2', 'line through last point with slope = -1/4', 'linear fit');
    p5fN.LineWidth = 3;
    p5fN2.LineWidth = 2;
    p5fN1.LineWidth = 2;
    p5fNh.LineWidth = 2;
    p5fNf.LineWidth = 2;
    p5fNfit.LineWidth = 2.5;
    p5fNfit.Color = [0.7, 0.7, 0.2];
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    grid on; 
    hold off;

fprintf('Template spacing for Tobs = %f hr: dfdot = %e * i^(%f)\n', Tobs_hr, 10^(templateFitCoeffs(2)), templateFitCoeffs(1));
end

%% Plot Maximum SNRg_i vs Tobs

if (PlotMaxSNRgvsTobs == 1)
    % Read in and process data from csv
    csvData = readmatrix('maxSNRgData.csv', 2, 0);
    Tobsvec = csvData(:, 1);
    Nsegvec = csvData(:, 2);
    imaxData = csvData(:, 3:2:end);
    SNRg_imaxData = csvData(:, 4:2:end);

    imax_avg = mean(imaxData, 2);
    imax_err = std(imaxData, 0, 2)/length(imaxData(1, :));
    SNRg_imax_avg = mean(SNRg_imaxData, 2);
    SNRg_imax_err = std(SNRg_imaxData, 0, 2)/length(SNRg_imaxData(1, :));

    figure(6+length(dfdot_sig)+N)
    subplot(1, 2, 1)
    e6fN1 = errorbar(Tobsvec, imax_avg, imax_err, '-k');
    hold on;
    p6fN1 = plot(Tobsvec, Nsegvec, '--b');
    title('i of Maximum SNRg_i vs Tobs');
    xlabel('Tobs (hr)');
    ylabel('i');
    legend('i of max SNRg', 'Nseg', 'Location', 'NorthWest');
    e6fN1.LineWidth = 2;
    grid on;
    ax = gca;
    ax.FontSize = 15;
    ax.LineWidth = 3;

    subplot(1, 2, 2)
    e6fN2 = errorbar(Tobsvec, SNRg_imax_avg, SNRg_imax_err, '-k');
    title('Maximum SNRg_i vs Tobs');
    xlabel('Tobs (hr)');
    ylabel('SNRg');
    e6fN2.LineWidth = 2;
    grid on;
    ax = gca;
    ax.FontSize = 15;
    ax.LineWidth = 3;
end
