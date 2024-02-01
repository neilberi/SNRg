%% Run TimeTemplates.m repeatedly

clear;
clc;
close all;

% Choose settings (0 = off, 1 = on)
ParsevalSNR = 0; % Nonsensical if on
NequalNseg = 0;

TS = 'f'; % Choose 'f' or 'fdot', but it should not matter for template times
TimeType = 'cpu'; % Choose 'cpu' or 'real'
OutputFile = 1;

%% Initialize and time template generation and SNRg calculation

filename = 'SNRgTemplateComputationTimes2.csv';
fout = fopen(filename, 'a');

Ntrials = 10;

fdotIn_vec = -5*10.^(-9:-5);
for it1 = 1:length(fdotIn_vec)
    if (fdotIn_vec(it1) == -5.e-9)
        TobsIn_vec = [8, 12, 16, 20, 24, 28, 32, 36, 40];
    else
        TobsIn_vec = [4, 16./3, 8, 12, 16, 20, 24, 28, 32];
    end

    for it2 = 1:length(TobsIn_vec)
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
        f_sig = 10.;
        
        % Signal frequency derivative (Hz/s)
        fdot_sig = fdotIn_vec(it1);
        
        % Length of observation (hr)
        Tobs_hr = TobsIn_vec(it2);
        Tobs = Tobs_hr * 3600.;
        
        % Coherence time (hr) - choose so that signal drifts 0.5 bin per coherence time
        if (abs(fdot_sig) > 0)
            Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
        else
            Tcoh_hr = 1./3600.*sqrt(0.5/abs(-5.e-6));
        end
        Tcoh = Tcoh_hr * 3600.;
        
        % Adjust coherence time to be an integer number of seconds
        Tcoh = floor(Tcoh);
        Tcoh_hr = Tcoh/3600.;
        
        % Adjust observation time to be exactly an integer number of coherence times to avoid binning headaches
        Nseg = floor(Tobs/Tcoh);
        Tobs = Nseg*Tcoh;
        Tobs_hr = Tobs / 3600.;
        fprintf('Observation time = %d sec (%f hr)\n',Tobs,Tobs_hr);
        fprintf('Coherence time per segment = %d sec (%f hr)\n',Tcoh,Tcoh_hr);
        
        % Calculate template spacings from model
        SNRgModelParams = readmatrix('SNRgModelParams.csv');
        TS_model = @(a, i, Tobs, fdot) 10^a(4)*Tobs^a(1)*i^a(2)*abs(fdot)^a(3);
        
        % Create vector of f/fdot values (for search)
        searchScale = 0;
        fStepSize = TS_model(SNRgModelParams(1, :), (Nseg+1)/2, Tobs_hr, fdot_sig);
        if (abs(fdot_sig) > 0)
            fdotStepSize = TS_model(SNRgModelParams(2, :), (Nseg+1)/2, Tobs_hr, fdot_sig);
        else
            fdotStepSize = TS_model(SNRgModelParams(2, :), (Nseg+1)/2, Tobs_hr, -5.e-6);
        end
        if (strcmp(TS, 'f'))
            fvec = (f_sig - searchScale*fStepSize):fStepSize:(f_sig + searchScale*fStepSize);
            fdotvec = fdot_sig;
            SNRg_iArray = zeros(Nseg, length(fvec));
        else
            fvec = f_sig;
            fdotvec = (fdot_sig - searchScale*fdotStepSize):fdotStepSize:(fdot_sig + searchScale*fdotStepSize);
            SNRg_iArray = zeros(Nseg, length(fdotvec));
        end
        
        % Set noise power mean from Parseval's Theorem
        noiseMeanPower = Tcoh*fsamp*noiseamp^2;
        
        % Define time series to hold raw data stream of signal plus noise
        deltat = 1./fsamp;
        t = 0:deltat:Tobs;
        Nsample = length(t);
        
        % Search band width (Hz)
        bandscale = 115; % Increase if fIndex is too high/low
        f_siglo = min(f_sig,f_sig+fdot_sig*Tobs);
        f_sighi = max(f_sig,f_sig+fdot_sig*Tobs);
        freqlo_approx = min(fvec)+min(fdotvec)*Tobs - (noiseOffset + nNoiseSample*(2*nBinSide+1))/Tcoh - bandscale/Tcoh;
        freqhi_approx = max(fvec) + max(fdotvec)*Tobs*(max(fdotvec) > 0) + (noiseOffset + nNoiseSample*(2*nBinSide+1))/Tcoh + bandscale/Tcoh;
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
        for seg = 1:Nseg
           fprintf('Generating segment %d\n',seg);
           indlo = (seg-1)*Nsample_coh + 1;
           indhi = indlo + Nsample_coh - 1;
           segment = data(indlo:indhi);
           segment_noNoise = signal(indlo:indhi);
           rawfft = fft(segment,Nsample_coh);
           spectrogram(1:nbandbin,seg) = abs(rawfft(indbandlo:indbandhi));
           rawSpectrogram(1:nbandbin,seg) = rawfft(indbandlo:indbandhi);
           [dummyVariable, fIndex_sig] = max(rawfft(indbandlo:indbandhi));
           if (seg == round(Nseg/5) || seg == round(Nseg/5*2) || seg == round(Nseg/5*3) || seg == round(Nseg/5*4) || seg == Nseg)
               fprintf('Signal fIndex of segment %d = %d; magnitude = %d; phase = %d\n', seg, fIndex_sig, abs(dummyVariable), angle(dummyVariable));
           end
        end

         % Plot spectrogram bar3(spectrogram)
        segarray = 1:Nseg;
        seghour = (segarray-1)*Tcoh/3600.;
        indarray = [indbandlo:indbandhi];
        freqplot = (indarray-indbandlo)*1.0/Tcoh + freqlo;

        % Create vector of indices for frequency trajectory in the middle of each time segment
        freqTrajIndex = round(((1:Nseg) - 0.5).*Nsample_coh) + 1;
        tMid = t(freqTrajIndex);
        
        % Calculate predicted nosie DS
        noiseMeanDS = noise_mean_DS(noiseMeanPower, 1, Nseg, Nseg);
        noiseStdDS = noise_std_DS(noiseMeanPower, 1, Nseg, Nseg);

        for it3 = 1:Ntrials
            f = @() GenerateTemplate2(fdot_sig, rawSpectrogram, f_sig, freqlo, freqhi, hamp, Tcoh, t, tMid, Nseg, nbandbin, indbandlo, indbandhi, Nsample_coh, nNoiseSample, nBinSide, noiseOffset, noiseMeanDS, noiseStdDS, ParsevalSNR);
            time = timeit(f);
            if (OutputFile == 1)
                fprintf(fout, '%e, %e, %e\n', fdot_sig, Tobs_hr, time);
            end
        end

    end
end

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