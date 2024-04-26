%% Script to test the effectiveness of interpolating between templates

clear;
clc;
close all;


CompS1S2 = 0;
CompS1S3 = 0;
CompS2S2tilde = 0;
CompS2S2tilde_sig = 1;
CompS2S2tilde_side = 0;
CompS1S3_sig = 0;
CompS2S2tildedots = 0;

Save = 0;

%% Initialize Templates

% Set number of template steps to interpolate
N_interp = 9;

% Choose interpolated templates ones to plot 
plotInds = 5; %[1, round(N_interp/2), N_interp];

% Set side bin number
nBinSide = 4;

% Set sampling frequency
f_samp = 32;

% Select initial f and fdot
f1 = 10;
%fdot1 = -5.e-6;

for fdot1 = -5.*10.^(-9:-5)
    if (fdot1 == -5*10^(-9))
        Tobs_hrvec = [8, 12, 16, 20, 24, 28, 32, 36, 40];
    else
        Tobs_hrvec = [4, 16./3, 8, 12, 16, 20, 24. 28, 32];
    end
for Tobs_hrSpace = Tobs_hrvec

% Select observation period and segmentation and based on first fdot
Tobs_hr = 40;
Tobs = 3600*Tobs_hr;
Tcoh = 1/sqrt(2*abs(fdot1));
Tcoh = round(Tcoh);
Nseg = floor(Tobs/Tcoh);
Tobs = Tcoh*Nseg;
Tobs_hr = Tobs/3600;

TobsSpace = 3600*Tobs_hrSpace;
NsegSpace = floor(TobsSpace/Tcoh);
TobsSpace = Tcoh*NsegSpace;
Tobs_hrSpace = TobsSpace/3600;

for g = 1:(2 - (fdot1==-5*10^(-9)) + 13*(fdot1==-5*10^(-5))):NsegSpace
% Select g that determines spacing
%g = 455;

% Load model parameters and calculate f and fdot for next template;
params = readmatrix('SNRgModelParams.csv');

TS_func = @(a, Tobs, g, fdot) 10.^( a(1).*log10(Tobs_hr) + a(2).*log10(g) + a(3).*log10(abs(fdot)) + a(4) );
fdot2 = zeros(1, N_interp);
%fdot2(1) = fdot1 + TS_func(params(2, :), Tobs_hr, g, fdot1);
fdot2(1) = fdot1 + TS_func(params(2, :), Tobs_hrSpace, g, fdot1);
for i = 2:N_interp
    %fdot2(i) = fdot2(i-1) + TS_func(params(2, :), Tobs_hr, g, fdot2(i-1));
    fdot2(i) = fdot2(i-1) + TS_func(params(2, :), Tobs_hrSpace, g, fdot2(i-1));
end

% Calculate band indices and frequencies
fdotvec = [fdot1, fdot2];

bandscale = 25; % Increase if fIndex is too high/low
f_siglo = min(f1,f1 + fdot1*Tobs);
f_sighi = max(f1,f1 + fdot1*Tobs);
freqlo_approx = f1 + min(fdotvec)*Tobs - bandscale/Tcoh;
freqhi_approx = f1 + bandscale/Tcoh;
freqhi = ceil(freqhi_approx * Tcoh)/Tcoh;
freqlo = floor(freqlo_approx * Tcoh)/Tcoh;
bandwidth = freqhi - freqlo;
indbandlo = floor(freqlo*Tcoh)+1;
indbandhi = floor(freqhi*Tcoh)+1;
nbandbin = indbandhi - indbandlo + 1;
segarray = 1:Nseg;
seghour = (segarray-1)*Tcoh/3600.;
indarray = indbandlo:indbandhi;
freqplot = (indarray-indbandlo)*1.0/Tcoh + freqlo;
fprintf('Low/high frequencies of band shown: %f-%f Hz\n',freqlo,freqhi);

%% Compute and Compare Templates

% Calculate SFT coefficients
S1 = GenerateTemplateCoeffs(f1, fdot1, Tobs, f_samp, Tcoh, indbandlo, indbandhi);
S2 = cell(1, N_interp);
for i = 1:N_interp
    S2{i} = GenerateTemplateCoeffs(f1, fdot2(i), Tobs, f_samp, Tcoh, indbandlo, indbandhi);
end


for i = plotInds
    if (~CompS1S2)
        break;
    end
    % Plot spectrograms
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, abs(S1).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, abs(S2{i}).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    % Plot magnitude ratios and phase differences
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, (abs(S2{i})./abs(S1)).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('%s/%s:\n%s = %e Hz/s\n%s = %e Hz/s\n', sprintf('$|S_{%d}|^2$', i+1), '$|S_1|^2$', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i), '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, phasediff(S2{i}, S1));
    s11.EdgeAlpha = 0.5;
    title(sprintf('Phase(%s)-Phase(%s):\n%s = %e Hz/s\n%s = %e Hz/s\n', sprintf('$S_{%d}$', i+1), '$S_1$', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i), '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
end

%% Generate third template and linearly interpolate magnitudes and phases between S3 and S1

% Calculate fdot3
%fdot3 = fdot2(end) + TS_func(params(2, :), Tobs_hr, g, fdot2(end));
fdot3 = fdot2(end) + TS_func(params(2, :), Tobs_hrSpace, g, fdot2(end));

% Calculate S3
S3 = GenerateTemplateCoeffs(f1, fdot3, Tobs, f_samp, Tcoh, indbandlo, indbandhi);

if (CompS1S3)
    % Plot spectrograms to compare S1 and S3
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, abs(S1).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, abs(S3).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', sprintf('%s%d}$', '$\dot{f}_{', N_interp+2), fdot3), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    % Plot magnitude ratios and phase differences
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, (abs(S3)./abs(S1)).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('%s/%s:\n%s = %e Hz/s\n%s = %e Hz/s\n', sprintf('$|S_{%d}|^2$', N_interp+2), '$|S_1|^2$', sprintf('%s%d}$', '$\dot{f}_{', N_interp+2), fdot3, '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, phasediff(S2{i}, S1));
    s11.EdgeAlpha = 0.5;
    title(sprintf('Phase(%s)-Phase(%s):\n%s = %e Hz/s\n%s = %e Hz/s\n', sprintf('$S_{%d}$', N_interp+2), '$S_1$', sprintf('%s%d}$', '$\dot{f}_{', N_interp+2), fdot3, '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
end


% Interpolate magnitudes and phases to approximate S2
S2_interpMag = cell(1, N_interp);
S2_interpPhase = cell(1, N_interp);
S2_interp = cell(1, N_interp);
for i = 1:N_interp
    S2_interpMag{i} = zeros(size(S2));
    S2_interpPhase{i} = zeros(size(S2));
    for t = 1:Nseg
        for k = 1:nbandbin
            S2_interpMag{i}(k, t) = interp1([fdot1, fdot3], [abs(S1(k, t)).^2, abs(S3(k, t)).^2], fdot2(i));
            %S2_interpPhase{i}(k, t) = interp1([fdot1, fdot3], [angle(S1(k, t)), angle(S3(k, t))], fdot2(i));
            S2_interpPhase{i}(k, t) = interp1Phase(fdot1, fdot3, angle(S1(k, t)), angle(S3(k, t)), fdot2(i));
        end
    end
    S2_interp{i} = sqrt(S2_interpMag{i}).*exp(sqrt(-1).*S2_interpPhase{i});
end

for i = plotInds
    if (~CompS2S2tilde)
        break;
    end
    % Plot spectrograms 
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, abs(S2{i}).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, abs(S2_interp{i}).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('Interpolated Template SFT Magnitudes:\n%s = %e Hz/s', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    % Plot magnitude ratios and phase differences
    figure;
    subplot(1, 2, 1);
    s11 = surf(seghour,freqplot, (abs(S2_interp{i})./abs(S2{i})).^2);
    s11.EdgeAlpha = 0.5;
    title(sprintf('%s/%s:\n%s = %e Hz/s\n', sprintf('%s%d}$', '$|\tilde{S_{', i+1), sprintf('$|S_{%d}|^2$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
    
    subplot(1, 2, 2);
    s11 = surf(seghour,freqplot, phasediff(S2_interp{i}, S2{i}));
    s11.EdgeAlpha = 0.5;
    title(sprintf('Phase(%s)-Phase(%s):\n%s = %e Hz/s\n', sprintf('%s%d}}$', '$\tilde{S_{',i+1), sprintf('$S_{%d}$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'LaTex');
    xlabel('Time (hr)', 'Interpreter', 'LaTex');
    ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
    view(0,90);
    ax = gca;
    ax.FontSize = 20;
    ax.LineWidth = 3;
    axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
    colorbar;
end

%% Do mag and phase plots as line plots along freqtraj and side bins

% Calculate frequency evolution of second template for the midpoint of each segment
tMid = ((1:Nseg) - 0.5)*Tcoh;
freqTraj2 = cell(1, N_interp);
for i = 1:N_interp
    freqTraj2{i} = f1 + fdot2(i).*tMid;
end
freqTraj1 = f1 + fdot1.*tMid;
freqTraj3 = f1 + fdot3.*tMid;

% Calculate frequency indices of signal
fIndex2 = cell(1, N_interp);
for i = 1:N_interp
    fIndex2{i} = round((freqTraj2{i}-freqlo)*Tcoh) + 1;
end
fIndex1 = round((freqTraj1-freqlo)*Tcoh) + 1;
fIndex3 = round((freqTraj1-freqlo)*Tcoh) + 1;

% Extract SFT coefficients of signal and side bins
S_sig1 = cell(1, N_interp);
S_sig2 = cell(1, N_interp);
S_sig3 = cell(1, N_interp);
S_side2 = cell(1, N_interp);
S_sig2_interp = cell(1, N_interp);
S_side2_interp = cell(1, N_interp);
for i = 1:N_interp
    S_sig1{i} = zeros(1, Nseg);
    S_sig2{i} = zeros(1, Nseg);
    S_sig3{i} = zeros(1, Nseg);
    S_side2{i} = zeros(2*nBinSide, Nseg);
    S_sig2_interp{i} = zeros(1, Nseg);
    S_side2_interp{i} = zeros(2*nBinSide, Nseg);
    for k = 1:Nseg
        S_sig1{i}(k) = S1(fIndex2{i}(k), k);
        S_sig2{i}(k) = S2{i}(fIndex2{i}(k), k);
        S_sig3{i}(k) = S3(fIndex2{i}(k), k);
        S_sig2_interp{i}(k) = S2_interp{i}(fIndex2{i}(k), k);
    
        for j = 1:nBinSide
            S_side2{i}(nBinSide-j+1, k) = S2{i}(fIndex2{i}(k)-j, k);
            S_side2{i}(nBinSide+j, k) = S2{i}(fIndex2{i}(k)+j, k);
            S_side2_interp{i}(nBinSide-j+1, k) = S2_interp{i}(fIndex2{i}(k)-j, k);
            S_side2_interp{i}(nBinSide+j, k) = S2_interp{i}(fIndex2{i}(k)+j, k);
        end
    end
end

% Calculate phase differences
segs = 1:Nseg;
dphi2_sig = cell(1, N_interp);
dphi2_side = cell(N_interp, 2*nBinSide);
for i = 1:N_interp
    dphi2_sig{i} = phasediff(S_sig2_interp{i}, S_sig2{i});
    for j = 1:2*nBinSide
        dphi2_side{i, j} = phasediff(S_side2_interp{i}(j, :), S_side2{i}(j, :));
    end
end

maxSeg = find(abs(phasediff(S_sig2_interp{(N_interp+1)/2}, S_sig2{(N_interp+1)/2}))>0.5, 1);
if (Save)
    fout = fopen(sprintf('MaxSegInterp%d.csv', (N_interp+1)/2), 'a');
    fprintf(fout, '%e, %e, %e, %e\n', fdot1, Tobs_hrSpace, g, maxSeg);
end

end
end
end


for i = plotInds
    if (~CompS2S2tilde_sig)
        break;
    end
    % Plot magnitude and phase differences of signal bins
    figure;
    subplot(2, 1, 1);
    p1 = plot(1:Nseg, (abs(S_sig2_interp{i})./abs(S_sig2{i})).^2, '-r');
    title(sprintf('Signal SFT: %s/%s\n%s = %e', sprintf('%s%d}}|^2$', '$|\tilde{S_{', i+1), sprintf('$|S_{%d}|^2$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
    p1.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
    
    subplot(2, 1, 2);
    p2 = plot(1:Nseg, phasediff(S_sig2_interp{i}, S_sig2{i}), '-b');
    title(sprintf('Signal SFT: Phase(%s) - Phase(%s)\n%s = %e', sprintf('%s%d}}$', '$\tilde{S_{',i+1), sprintf('$S_{%d}$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Phase Difference', 'Interpreter', 'Latex');
    p2.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
end

for i = plotInds
    if (~CompS1S3_sig)
        break;
    end
    % Plot magnitude and phase differences of signal bins for S1 and S3
    figure;
    subplot(2, 1, 1);
    p1 = plot(1:Nseg, (abs(S_sig3{i})./abs(S_sig1{i})).^2, '-r');
    title(sprintf('$S_%d$ Signal SFT: %s/%s\n%s = %e\n%s = %e\n%s = %e', i+1, sprintf('%s%d}}|^2$', '$|\tilde{S_{', N_interp+2), sprintf('$|S_{%d}|^2$', 1), sprintf('%s%d}$', '$\dot{f}_{', N_interp+2), fdot3, '$\dot{f}_1$', fdot1), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
    p1.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
    
    subplot(2, 1, 2);
    p2 = plot(1:Nseg, phasediff(S_sig3{i}, S_sig1{i}), '-b');
    title(sprintf('$S_%d$ Signal SFT: Phase(%s) - Phase(%s)\n%s = %e\n%s = %e', i+1, sprintf('%s%d}}$', '$\tilde{S_{',N_interp+2), sprintf('$S_{%d}$', 1), sprintf('%s%d}$', '$\dot{f}_{', N_interp+2), fdot3, '$\dot{f}_1$', fdot1), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Phase Difference', 'Interpreter', 'Latex');
    p2.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
end

for i = plotInds
    if (~CompS2S2tilde_side)
        break;
    end
    % Plot magnitude and phase differences of side bins
    figure;
    sgtitle(sprintf('Side Bin SFT:\n%s = %e', sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'Latex');
    for j = 1:2*nBinSide
        subplot(2, 2*nBinSide, j);
        p1 = plot(1:Nseg, (abs(S_side2_interp{i}(j, :))./abs(S_side2{i}(j, :))).^2, '-r');
        title(sprintf('Side Bin %d: %s/%s', j, sprintf('%s%d}}|^2$', '$|\tilde{S_{', i+1), sprintf('$|S_{%d}|^2$', i+1)), 'Interpreter', 'Latex');
        xlabel('Time Segment', 'Interpreter', 'Latex');
        ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
        p1.LineWidth = 2;
        ax = gca;
        ax.LineWidth = 3;
        ax.FontSize = 8;
        grid on;
        
        subplot(2, 2*nBinSide, j+2*nBinSide);
        p2 = plot(1:Nseg, phasediff(S_side2_interp{i}(j, :), S_side2{i}(j, :)), '-b');
        title(sprintf('Side Bin %d: Phase(%s) - Phase(%s)', j, sprintf('%s%d}}$', '$\tilde{S_{',i+1), sprintf('$S_{%d}$', i+1)), 'Interpreter', 'Latex');
        xlabel('Time Segment', 'Interpreter', 'Latex');
        ylabel('Phase Difference', 'Interpreter', 'Latex');
        p2.LineWidth = 2;
        ax = gca;
        ax.LineWidth = 3;
        ax.FontSize = 8;
        grid on;
    end
end


%% Calculate Signal Dot Products for S2 and ~S2

% Create noise
nbandbin = indbandhi-indbandlo+1;
Nsample_coh = f_samp*Tcoh;
Nseg = floor(Tobs/Tcoh);
N_sample = length(0:(1/f_samp):Tobs);
noiseamp = 12;
noise = noiseamp*random('norm',0.,1.,1., N_sample);
n = (1+sqrt(-1)).*zeros(nbandbin,Nseg);
for seg = 1:Nseg
    indlo = (seg-1)*Nsample_coh + 1;
    indhi = indlo + Nsample_coh - 1;
    nSegment = noise(indlo:indhi);
    nRawfft = fft(nSegment, Nsample_coh);
    n(1:nbandbin,seg) = nRawfft(indbandlo:indbandhi);
end

Sdots2 = cell(N_interp, 1);
Sdots2_interp = cell(N_interp, 1);
for i = 1:N_interp
    Sdots2{i} = zeros(1, Nseg);
    Sdots2_interp{i} = zeros(1, Nseg);
    for tIndex = 1:Nseg
        weight = 1/sqrt(sum(abs( S2{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ).^2 )));
        Sdots2{i}(tIndex) = sum(conj( S2{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ) ).*( S2{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex )) + n( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ) )/weight;

        weight_interp = 1/sqrt(sum(abs( S2_interp{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ).^2 )));
        Sdots2_interp{i}(tIndex) = sum(conj( S2_interp{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ) ).*( S2{i}( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ))  + n( (fIndex2{i}(tIndex)-nBinSide):(fIndex2{i}(tIndex)+nBinSide), tIndex ) )/weight_interp;
    end
end

for i = plotInds
    if (~CompS2S2tildedots)
        break;
    end
    % Plot magnitude and phase differences of signal bins
    figure;
    subplot(2, 1, 1);
    p1 = plot(1:Nseg, (abs(Sdots2_interp{i})./abs(Sdots2{i})).^2, '-r');
    title(sprintf('Signal Dot Products: %s/%s\n%s = %e', sprintf('%s%d}}|^2$', '$|\tilde{S_{', i+1), sprintf('$|S_{%d}|^2$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
    p1.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
    
    subplot(2, 1, 2);
    p2 = plot(1:Nseg, phasediff(Sdots2_interp{i}, Sdots2{i}), '-b');
    title(sprintf('Signal Dot Products: Phase(%s) - Phase(%s)\n%s = %e', sprintf('%s%d}}$', '$\tilde{S_{',i+1), sprintf('$S_{%d}$', i+1), sprintf('%s%d}$', '$\dot{f}_{', i+1), fdot2(i)), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Phase Difference', 'Interpreter', 'Latex');
    p2.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 18;
    grid on;
end
