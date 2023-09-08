% Simplified CalculateSNRg.m with timers to help calculate computational
% cost of search

clear;
clc;
close all;
profile clear

% Choose settings (0 = off, 1 = on)
ParsevalSNR = 1;
NequalNseg = 0;
Skip2TimeFields = 0;

TS = 'f'; % Choose 'f' or 'fdot'
TimeType = 'real'; % Choose 'cpu' or 'real'

%% Generate Data Spectrogram

% SNR1 and {SNRg_i | M<=i<=N} will be calculated
M = 3;
N = 3;

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
Tobs_hr = 8.;
Tobs = Tobs_hr * 3600.;

% Coherence time (hr) - choose so that signal drifts 0.5 bin per coherence time
Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
Tcoh = Tcoh_hr * 3600.;

% Adjust coherence time to be an integer number of seconds
Tcoh = floor(Tcoh);
Tcoh_hr = Tcoh/3600.;

% DO NOT TOUCH/CARE
searchScale = 0;

fStepSize = (0.1)/Tcoh*(Tcoh < Tobs) + 0.3/sqrt(0.5/(1.e-6))*(Tcoh >= Tobs);
fvec = (f_sig-searchScale*fStepSize) : fStepSize : (f_sig+searchScale*fStepSize);

fdotStepSize = (7.8e-6)*abs(fdot_sig)*(fdot_sig ~= 0) + 1.e-7*(fdot_sig == 0);
fdotvec = (fdot_sig-searchScale*fdotStepSize) : fdotStepSize : (fdot_sig+searchScale*fdotStepSize);

[F, Fdot] = meshgrid(fvec, fdotvec);

% Adjust observation time to be exactly an integer number of coherence timestruct to avoid binning headaches
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
freqlo_approx = min(fvec)+min(fdotvec)*Tobs - (noiseOffset + nNoiseSample + nBinSide)/Tcoh - bandscale/Tcoh;
freqhi_approx = max(fvec)+max(fdotvec)*Tobs*(max(fdotvec) > 0) + (noiseOffset + nNoiseSample + nBinSide)/Tcoh + bandscale/Tcoh;
freqhi = ceil(freqhi_approx * Tcoh)/Tcoh;
freqlo = floor(freqlo_approx * Tcoh)/Tcoh;
bandwidth = freqhi - freqlo;
fprintf('Low/high frequencies of band shown: %f-%f Hz\n',freqlo,freqhi);

% Number of bins signal drifts per coherence time
Nbin_drift = abs(fdot_sig) * Tcoh * Tcoh;
Nbin_drift_total = abs(fdot_sig) * Tobs * Tcoh;
fprintf('Number of bins signal drifts per coherence time = %8.2f\n',Nbin_drift);
fprintf('Total number of bins signal drifts = %8.2f\n',Nbin_drift_total);

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
s1.EdgeAlpha = 0.5;
title('Raw spectrogram');
xlabel('Time (hours)');
ylabel('Frequency (Hz)');
view(0,90);
ax = gca;
ax.FontSize = 18;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);

%% Generate Templates and Calculate SNRs
% Calculate template spacings from model
SNRgModelParams = readmatrix('SNRgModelParams.csv');
TS_model = @(a, i, Tobs_hr, fdot) 10^a(4)*Tobs^a(1)*i^a(2)*abs(fdot)^a(3);

% Create vector of f/fdot values (for search)
searchScale = 1;
fStepSize = TS_model(SNRgModelParams(1, :), (N+M)/2, Tobs_hr, fdot_sig);
fdotStepSize = TS_model(SNRgModelParams(2, :), (N+M)/2, Tobs_hr, fdot_sig);
if (strcmp(TS, 'f'))
    fvec = (f_sig - searchScale*fStepSize):fStepSize:(f_sig + searchScale*fStepSize);
    fdotvec = fdot_sig;
    SNR1Array = zeros(1, length(fvec));
    SNRg_iArray = zeros(N-M+1, length(fvec));
else
    fvec = f_sig;
    fdotvec = (fdot_sig - searchScale*fdotStepSize):fdotStepSize:(fdot_sig + searchScale*fdotStepSize);
    SNR1Array = zeros(1, length(fdotvec));
    SNRg_iArray = zeros(N-M+1, length(fdotvec));
end

% Create vector of indices for frequency trajectory in the middle of each time segment
freqTrajIndex = ((1:Nseg) - 0.5).*Nsample_coh + 1;
tMid = t(freqTrajIndex);

% Create arrays to store profiler info structs
profiles = cell(length(fvec), length(fdotvec));

for r = 1:length(fvec)
    for j = 1:length(fdotvec)
        
        if (strcmp(TimeType, 'cpu'))
            profile on -timer 'cpu';
        else
            profile on -timer 'real';
        end
        
        if (strcmp(TS, 'f'))
            fprintf('Calulating SNRg_is for df = %e Hz\n', fvec(r)-f_sig);
        else
            fprintf('Calulating SNRg_is for dfdot = %e Hz/s\n', fdotvec(j)-fdot_sig);
        end
   
        % Calculate expected phase/frequency trajectories and signal
        freqTraj = fvec(r) + fdotvec(j).*tMid; % Nseg add, Nseg mult
        phaseTraj = 2*pi.*(fvec(r).*t + 0.5*fdotvec(j).*t.^2); % Nsample add, 4*Nsample mult
        searchSignal = hamp.*sin(phaseTraj); % Nsample sin
    
        % Create template spectrogram
        searchSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
        for seg = 1:Nseg
            indlo = (seg-1)*Nsample_coh + 1; % 2*Nseg add, Nseg mult
            indhi = indlo + Nsample_coh - 1; % 2*Nseg add
            searchSegment = searchSignal(indlo:indhi); 
            searchRawfft = fft(searchSegment, Nsample_coh); % Nseg ffts of Nsample_coh length data sets
            searchSpectrogram(1:nbandbin,seg) = searchRawfft(indbandlo:indbandhi)*fftNormalization; % Nseg*nbandbin mult
        end
    
        % Calculate signal frequency indices
        fIndex = round((freqTraj-freqlo)*Tcoh) + 1; % 2*Nseg add, Nseg mult
        
        % Calculate background noise frequency indices
        if (ParsevalSNR == 0)
            fIndex_noise = zeros(nNoiseSample, Nseg);
            for k = 1:nNoiseSample
                fIndex_noise(k, :) = fIndex + 10*k + noiseOffset;
            end
        end
    
        % Calculate weights array
        weight = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
        for tIndex = 1:Nseg
            weight(:, tIndex) = searchSpectrogram(:, tIndex).*abs(searchSpectrogram(:, tIndex)).^2/sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2); % 2*nBinSide*Nseg + 2*Nseg add, 3*Nseg*nbandbin + 2*(2*nBinSide+1)*Nseg mult
        end
    
        % Calculate detection statistics
        signalDS1 = 0;
        signalDSg_i = zeros(N-M+1, 1);
        grouper_sig = (1+sqrt(-1)).*zeros(N-M+1, 1);
        if (ParsevalSNR == 0)
            noiseDS1 = zeros(1, nNoiseSample);
            noiseDSg_i = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);
            grouper_noise = (1+sqrt(-1)).*zeros(N-M+1, nNoiseSample);
        end
    
        for tIndex = 1:Nseg
            % Signal detection statistics
            grouper_sig = grouper_sig + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)); % 5*Nseg add, (2*nBinSide+1)*Nseg mult
            for i = M:N  % ceil(Nseg/i) add, ceil(Nseg/i) mult
                if (mod(tIndex, i) == 0)
                    signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + abs(grouper_sig(i-M+1))^2;
                    grouper_sig(i-M+1) = 0;
                elseif (tIndex == Nseg)
                    signalDSg_i(i-M+1) = signalDSg_i(i-M+1) + abs(grouper_sig(i-M+1))^2;
                end
            end
    
            % Noise detection statistics
            if (ParsevalSNR == 0)
                for k = 1:nNoiseSample
                    grouper_noise(:, k) = grouper_noise(:, k) + sum(conj(weight((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex_noise(k, tIndex)-nBinSide):(fIndex_noise(k, tIndex)+nBinSide), tIndex));
                end
                for i = M:N
                    if (mod(tIndex, i) == 0)
                        noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + abs(grouper_noise(i-M+1, :)).^2;
                        grouper_noise(i-M+1, :) = 0;
                    elseif (tIndex == Nseg)
                        noiseDSg_i(i-M+1, :) = noiseDSg_i(i-M+1, :) + abs(grouper_noise(i-M+1, :)).^2;
                    end
                end
            end
        end
        
        % Calculate SNRs
        for i = M:N 
            if (ParsevalSNR == 0)
                SNRg_iArray(i-M+1, j*r) = abs(signalDSg_i(i-M+1) - mean(noiseDSg_i(i-M+1, :)))/std(noiseDSg_i(i-M+1, :));
            elseif (ParsevalSNR == 1)
                SNRg_iArray(i-M+1, j*r) = abs(signalDSg_i(i-M+1) - noiseMean)/noiseSTD/sqrt(i);
            end
        end
        
        profiles{r, j} = profile('info');
    end
end

%% Calculate and Organize Times

% Create structures to store times
timestruct.TimeDomainSignal = 0;  
timestruct.SegmentAndFFT = 0;
timestruct.CalcfIndex = 0;
timestruct.CalcWeights = 0;
timestruct.CalcDS = 0;
timestruct.CalcSNR = 0;
times = cell(length(fvec), length(fdotvec));

% Fill structs with times from profiler
for k = 1:numel(times)
    times{k} = timestruct;
    execLines = profiles{k}.FunctionTable.ExecutedLines;
    
    times{k}.TimeDomainSignal = sum(execLines(logical((209 <= execLines(:, 1)).*(execLines(:, 1) <= 211)), 3));
    times{k}.SegmentAndFFT = sum(execLines(logical((214 <= execLines(:, 1)).*(execLines(:, 1) <= 221)), 3));
    times{k}.CalcfIndex = sum(execLines(logical((224 <= execLines(:, 1)).*(execLines(:, 1) <= 232)), 3));
    times{k}.CalcWeights = sum(execLines(logical((235 <= execLines(:, 1)).*(execLines(:, 1) <= 238)), 3));
    times{k}.CalcDS = sum(execLines(logical((241 <= execLines(:, 1)).*(execLines(:, 1) <= 276)), 3));
    times{k}.CalcSNR = sum(execLines(logical((279 <= execLines(:, 1)).*(execLines(:, 1) <= 285)), 3));
end

% Plot times for each template and for all templates
fn = fieldnames(timestruct);
if (Skip2TimeFields == 1)
    fn = {fn{3}, fn{4}, fn{5}, fn{6}};
end
bars = zeros(numel(times), numel(fn));
legends = cell(1, numel(times));
cats = categorical(fn);
cats = reordercats(cats, fn);
for k = 1:numel(times)
    if (strcmp(TS, 'f'))
        legends{k} = sprintf('f = %f Hz', fvec(k));
    else
        legends{k} = sprintf('fdot = %e Hz/s', fdotvec(k));
    end
    for j = 1:length(fn)
        bars(k, j) = times{k}.(fn{j});
    end

    figure;
    bk = bar(cats, bars(k, :));
    title(['Run Time for Various Template Calculations for ', legends{k}]);
    ylabel([TimeType, ' time (s)']);
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 15;
    grid on;
end

figure;
b1 = bar(cats, bars, 'stacked');
title('Run Time for Various Template Calculations');
ylabel([TimeType, ' time (s)']);
legend(legends);
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 15;
grid on;


