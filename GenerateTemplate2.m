function [SNRg] = GenerateTemplate(fdot_sig, g, rawSpectrogram, f_sig, freqlo, freqhi, hamp, Tcoh, t, tMid, Nseg, nbandbin, indbandlo, indbandhi, Nsample_coh, nNoiseSample, nBinSide, noiseOffset, noiseMeanDS, noiseStdDS, ParsevalSNR)
%Generates template and calculates SNRg
fprintf('Calulating SNRg_is for df = %e Hz\n', 0);

% Calculate expected phase/frequency trajectories and signal
freqTraj = f_sig + fdot_sig.*tMid; 
phaseTraj = 2*pi.*(f_sig.*t + 0.5*fdot_sig.*t.^2);
searchSignal = hamp.*sin(phaseTraj);

% Create template spectrogram
searchSpectrogram = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
for seg = 1:Nseg
    indlo = (seg-1)*Nsample_coh + 1;
    indhi = indlo + Nsample_coh - 1;
    searchSegment = searchSignal(indlo:indhi);
    searchRawfft = fft(searchSegment, Nsample_coh);
    searchSpectrogram(1:nbandbin,seg) = searchRawfft(indbandlo:indbandhi);
end

% Calculate signal frequency indices
fIndex = round((freqTraj-freqlo)*Tcoh) + 1;

% Calculate background noise frequency indices
if (ParsevalSNR == 0)
    fIndex_noise = zeros(nNoiseSample, Nseg);
    for k = 1:nNoiseSample
        fIndex_noise(k, :) = fIndex + (2*nBinSide+2)*k + noiseOffset;
    end
end

% Calculate weightedFFTs array
weight = zeros(1, Nseg);
weightedFFT = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
for tIndex = 1:Nseg
    weight(tIndex) = 1/sqrt(sum(abs(searchSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).^2));
    weightedFFT(:, tIndex) = searchSpectrogram(:, tIndex).*weight(:, tIndex);
end

% Calculate detection statistics
signalDSg_i = zeros(1, 1);
grouper_sig = (1+sqrt(-1)).*zeros(1, 1);
if (ParsevalSNR == 0)
    noiseDSg_i = (1+sqrt(-1)).*zeros(1, nNoiseSample);
    grouper_noise = (1+sqrt(-1)).*zeros(1, nNoiseSample);
end

for tIndex = 1:Nseg
    % Signal detection statistics
    grouper_sig = grouper_sig + sum(conj(weightedFFT((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex));
    for i = g:g

        contribution = abs(grouper_sig(1))^2;

        if (mod(tIndex, i) == 0)
            signalDSg_i(1) = signalDSg_i(1) + contribution;
            grouper_sig(1) = 0;
        elseif (tIndex == Nseg)
            signalDSg_i(1) = signalDSg_i(1) + contribution;
        end
    end

    % Noise detection statistics
    if (ParsevalSNR == 0)
        for k = 1:nNoiseSample
            grouper_noise(:, k) = grouper_noise(:, k) + sum(conj(weightedFFT((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex_noise(k, tIndex)-nBinSide):(fIndex_noise(k, tIndex)+nBinSide), tIndex));
        end
        for i = g:g

            contribution = abs(grouper_noise(1, :)).^2;

            if (mod(tIndex, i) == 0)
                noiseDSg_i(1, :) = noiseDSg_i(1, :) + contribution;
                grouper_noise(1, :) = 0;
            elseif (tIndex == Nseg)
                noiseDSg_i(1, :) = noiseDSg_i(1, :) + contribution;
            end
        end
    end
end


for i = g:g
    if (ParsevalSNR == 0)
        SNRg = abs(signalDSg_i(1) - mean(noiseDSg_i(1, :)))/std(noiseDSg_i(1, :));
    elseif (ParsevalSNR == 1)
        SNRg = abs(signalDSg_i(1) - noiseMeanDS(g))/noiseStdDS(g);
    end
end
end