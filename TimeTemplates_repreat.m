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

filename = 'SNRgTemplateComputationTimes.csv';
fout = fopen(filename, 'a');

Ntrials = 5;

fdotIn_vec = -5*10.^(-9:-5 );
for it1 = 1:length(fdotIn_vec)
    fdotIn = fdotIn_vec(it1);
    if (fdotIn == -5.e-9)
        TobsIn_vec = [8, 12, 16, 20, 24, 28, 32, 36, 40];
    else
        TobsIn_vec = [4, 16./3, 8, 12, 16, 20, 24, 28, 32];
    end

    for it2 = 1:length(TobsIn_vec)
        TobsIn = TobsIn_vec(it2);
        
        Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdotIn));
        Tcoh = Tcoh_hr * 3600.;
        Tcoh = floor(Tcoh);
        Nseg = floor(TobsIn*3600/Tcoh);
        
        for it3 = 1:Nseg
            iIn = it3;
            for trial = 1:Ntrials
    
                % SNR1 and {SNRg_i | M<=i<=N} will be calculated
                M = 1;
                N = 1;
                if (OutputFile == 1)
                    M = iIn;
                    N = iIn;
                end
                
                % Set information on noise samples
                nNoiseSample = 10;
                noiseOffset = 50;
                nBinSide = 4;
                
                % Sampling frequency
                fsamp = 128.;
                
                % Noise ASD ( / sqrt(Hz) )
                hnoise = 12;
                noiseamp = hnoise;
                
                % Signal initial frequency (Hz)
                f_sig = 10.;
                
                % Signal frequency derivative (Hz/s)
                fdot_sig = fdotIn;
                
                % Length of observation (hr)
                Tobs_hr = TobsIn;
                Tobs = Tobs_hr * 3600.;
                
                % Coherence time (hr) - choose so that signal drifts 0.5 bin per coherence time
                Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
                Tcoh = Tcoh_hr * 3600.;
                
                % Adjust coherence time to be an integer number of seconds
                Tcoh = floor(Tcoh);
                Tcoh_hr = Tcoh/3600.;
                
                % Adjust observation time to be exactly an integer number of coherence timestruct to avoid binning headaches
                Nseg = floor(Tobs/Tcoh);
                if (NequalNseg == 1)
                    N = Nseg;
                end
                
                Tobs = Nseg*Tcoh;
                Tobs_hr = Tobs / 3600.;
                fprintf('Observation time = %d sec (%f hr)\n',Tobs,Tobs_hr);
                fprintf('Coherence time per segment = %d sec (%f hr)\n',Tcoh,Tcoh_hr);
                assert(M<=N, 'Error: M>N');
                
                % Calculate template spacings from model
                SNRgModelParams = readmatrix('SNRgModelParams.csv');
                TS_model = @(a, i, Tobs_hr, fdot) 10^a(4)*Tobs^a(1)*i^a(2)*abs(fdot)^a(3);
                
                % Create vector of f/fdot values (for search)
                searchScale = 0;
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
                
                % Set fft normalizations
                % Set noise power mean and standard deviation from Parseval's Theorem
                fftNormalization = (Tcoh*fsamp*noiseamp^2)^(-0.5);
                if (ParsevalSNR == 1)
                    noiseMean = Nseg*(2*nBinSide+1)^2;
                    noiseSTD = Nseg*(2*nBinSide+1)^2;
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
                hamp = 1.;
                
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
                
                if (OutputFile == 0)
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
                            searchSpectrogram(1:nbandbin,seg) = searchRawfft(indbandlo:indbandhi); %*fftNormalization; % Nseg*nbandbin mult
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
                        if (ParsevalSNR == 1)
                            noiseMeanPower = Tcoh*fsamp*noiseamp^2;
                            noiseMeanDS = noise_mean_DS(noiseMeanPower, M, N, Nseg);
                            noiseStdDS = noise_std_DS(noiseMeanPower, M, N, Nseg);
                        end

                        for i = M:N 
                            if (ParsevalSNR == 0)
                                SNRg_iArray(i-M+1, j*r) = abs(signalDSg_i(i-M+1) - mean(noiseDSg_i(i-M+1, :)))/std(noiseDSg_i(i-M+1, :));
                            elseif (ParsevalSNR == 1)
                                SNRg_iArray(i-M+1, j*r) = abs(signalDSg_i(i-M+1) - noiseMeanDS(i-M+1))/noiseStdDS(i-M+1);
                            end
                        end
                        
                        profiles{r, j} = profile('info');
                    end
                end
                
                
                % Fill structs with times from profiler
                execLines = profiles{1}.FunctionTable.ExecutedLines;
                
               
                if (OutputFile == 1)
                    fprintf(fout, '%e, %e, %e, %e\n', Tobs_hr, N, fdot_sig, sum(execLines(:, 3)));
                end
            end

        end
    end
end

%% Functions

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