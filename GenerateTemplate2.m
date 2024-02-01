function [pow2s, SNRg] = GenerateTemplate2(fdot_sig, rawSpectrogram, f_sig, freqlo, freqhi, hamp, Tcoh, t, tMid, Nseg, nbandbin, indbandlo, indbandhi, Nsample_coh, nNoiseSample, nBinSide, noiseOffset, noiseMeanDS, noiseStdDS, ParsevalSNR)

pow2s = 2.^(1:floor(log2(Nseg)));
N2s = length(pow2s);

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
signalDSg = zeros(N2s, 1);
grouper_sig = (1+sqrt(-1)).*zeros(N2s, 1);
if (ParsevalSNR == 0)
    noiseDSg = (1+sqrt(-1)).*zeros(N2s, nNoiseSample);
    grouper_noise = (1+sqrt(-1)).*zeros(N2s, nNoiseSample);
end

for tIndex = 1:Nseg
    % Signal detection statistics
    grouper_sig = grouper_sig + sum(conj(weightedFFT((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex));
    for g = pow2s

        contribution = abs(grouper_sig(log2(g)))^2;

        if (mod(tIndex, g) == 0)
            signalDSg(log2(g)) = signalDSg(log2(g)) + contribution;
            grouper_sig(1) = 0;
        elseif (tIndex == Nseg)
            signalDSg(log2(g)) = signalDSg(log2(g)) + contribution;
        end
    end

    % Noise detection statistics
    if (ParsevalSNR == 0)
        for k = 1:nNoiseSample
            grouper_noise(:, k) = grouper_noise(:, k) + sum(conj(weightedFFT((fIndex(tIndex)-nBinSide):(fIndex(tIndex)+nBinSide), tIndex)).*rawSpectrogram((fIndex_noise(k, tIndex)-nBinSide):(fIndex_noise(k, tIndex)+nBinSide), tIndex));
        end
        for g = pow2s

            contribution = abs(grouper_noise(log2(g), :)).^2;

            if (mod(tIndex, g) == 0)
                noiseDSg(log2(g), :) = noiseDSg(log2(g), :) + contribution;
                grouper_noise(1, :) = 0;
            elseif (tIndex == Nseg)
                noiseDSg(log2(g), :) = noiseDSg(log2(g), :) + contribution;
            end
        end
    end
end

SNRg = zeros(N2s, 1);
for g = pow2s
    if (ParsevalSNR == 0)
        SNRg(log2(g)) = abs(signalDSg(log2(g)) - mean(noiseDSg(log2(g), :)))/std(noiseDSg(log2(g), :));
    elseif (ParsevalSNR == 1)
        SNRg(log2(g)) = abs(signalDSg(log2(g)) - noiseMeanDS(g))/noiseStdDS(g);
    end
end
end