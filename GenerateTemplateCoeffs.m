function [searchSpectrogram] = GenerateTemplateCoeffs(f, fdot, Tobs, f_samp, Tcoh, indbandlo, indbandhi, hamp)

if(~exist('hamp', 'var'))
    hamp = 1;
end

nbandbin = indbandhi-indbandlo+1;
Nsample_coh = f_samp*Tcoh;
Nseg = floor(Tobs/Tcoh);

t = 0:(1/f_samp):Tobs;

phaseTraj = 2*pi.*(f.*t + 0.5*fdot.*t.^2);
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

end