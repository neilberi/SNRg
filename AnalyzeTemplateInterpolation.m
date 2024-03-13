%% Script to test the effectiveness of interpolating between templates

clear;
clc;
close all;

%% Initialize Templates

% Set side bin number
nBinSide = 4;

% Set sampling frequency and observation period
f_samp = 32;

Tobs_hr = 24;
Tobs = 3600*Tobs_hr;

% Select g that determines spacing
g = 1;

% Select initial f and fdot
f1 = 10;
fdot1 = -5.e-6;

% Load model parameters and calculate f and fdot for next template;
params = readmatrix('SNRgModelParams.csv');

TS_func = @(a, Tobs, g, fdot) 10.^( a(1).*log10(Tobs) + a(2).*log10(g) + a(3).*log10(abs(fdot)) + a(4) );
fdot2 = fdot1 + TS_func(params(2, :), Tobs, g, fdot1);

% Select segmentation and based on first fdot
Tcoh = 1/sqrt(2*abs(fdot1));
Tcoh = round(Tcoh);
Nseg = floor(Tobs/Tcoh);
Tobs = Tcoh*Nseg;

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
S2 = GenerateTemplateCoeffs(f1, fdot2, Tobs, f_samp, Tcoh, indbandlo, indbandhi);

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
s11 = surf(seghour,freqplot, abs(S2).^2);
s11.EdgeAlpha = 0.5;
title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', '$\dot{f}_2$', fdot2), 'Interpreter', 'LaTex');
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
s11 = surf(seghour,freqplot, (abs(S2)./abs(S1)).^2);
s11.EdgeAlpha = 0.5;
title(sprintf('Magnitude(%s)/Magnitude(%s):\n%s = %e Hz/s\n%s = %e Hz/s\n', '$S_2$', '$S_1$', '$\dot{f}_2$', fdot2, '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
xlabel('Time (hr)', 'Interpreter', 'LaTex');
ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
view(0,90);
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
colorbar;

subplot(1, 2, 2);
s11 = surf(seghour,freqplot, phasediff(S2, S1));
s11.EdgeAlpha = 0.5;
title(sprintf('Phase(%s)-Phase(%s):\n%s = %e Hz/s\n%s = %e Hz/s\n', '$S_2$', '$S_1$', '$\dot{f}_2$', fdot2, '$\dot{f}_1$', fdot1), 'Interpreter', 'LaTex');
xlabel('Time (hr)', 'Interpreter', 'LaTex');
ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
view(0,90);
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
colorbar;

%% Generate third template and linearly interpolate magnitudes and phases between S3 and S1

% Calculate fdot3
fdot3 = fdot2 + TS_func(params(2, :), Tobs, g, fdot2);

% Calculate S3
S3 = GenerateTemplateCoeffs(f1, fdot3, Tobs, f_samp, Tcoh, indbandlo, indbandhi);

% Interpolate magnitudes and phases to approximate S2
S2_interpMag = zeros(size(S2));
S2_interpPhase = zeros(size(S2));
for t = 1:Nseg
    for k = 1:nbandbin
        S2_interpMag(k, t) = interp1([fdot1, fdot3], [abs(S1(k, t)).^2, abs(S3(k, t)).^2], fdot2);
        S2_interpPhase(k, t) = interp1([fdot1, fdot3], [phase(S1(k, t)), phase(S3(k, t))], fdot2);
    end
end
S2_interp = sqrt(S2_interpMag).*exp(sqrt(-1).*S2_interpPhase);

% Plot spectrograms 
figure;
subplot(1, 2, 1);
s11 = surf(seghour,freqplot, abs(S2).^2);
s11.EdgeAlpha = 0.5;
title(sprintf('Template SFT Magnitudes:\n%s = %e Hz/s', '$\dot{f}_2$', fdot2), 'Interpreter', 'LaTex');
xlabel('Time (hr)', 'Interpreter', 'LaTex');
ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
view(0,90);
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
colorbar;

subplot(1, 2, 2);
s11 = surf(seghour,freqplot, abs(S2_interp).^2);
s11.EdgeAlpha = 0.5;
title(sprintf('Interpolated Template SFT Magnitudes:\n%s = %e Hz/s', '$\dot{f}_2$', fdot2), 'Interpreter', 'LaTex');
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
s11 = surf(seghour,freqplot, (abs(S2_interp)./abs(S2)).^2);
s11.EdgeAlpha = 0.5;
title(sprintf('Magnitude(%s)/Magnitude(%s):\n%s = %e Hz/s\n', '$\tilde{S_2}$', '$S_2$', '$\dot{f}_2$', fdot2), 'Interpreter', 'LaTex');
xlabel('Time (hr)', 'Interpreter', 'LaTex');
ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
view(0,90);
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
colorbar;

subplot(1, 2, 2);
s11 = surf(seghour,freqplot, phasediff(S2_interp, S2));
s11.EdgeAlpha = 0.5;
title(sprintf('Phase(%s)-Phase(%s):\n%s = %e Hz/s\n', '$\tilde{S_2}$', '$S_2$', '$\dot{f}_2$', fdot2), 'Interpreter', 'LaTex');
xlabel('Time (hr)', 'Interpreter', 'LaTex');
ylabel('Frequency (Hz)', 'Interpreter', 'LaTex');
view(0,90);
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
axis([min(seghour) max(seghour) min(freqplot) max(freqplot)]);
colorbar;

%% Do mag and phase plots as line plots along freqtraj and side bins

% Set maximum dphi for 3rd plot
phimax = 0.01;

% Calculate frequency evolution of second template for the midpoint of each segment
tMid = ((1:Nseg) - 0.5)*Tcoh;
freqTraj2 = f1 + fdot2.*tMid;

% Calculate frequency indices of signal
fIndex2 = round((freqTraj2-freqlo)*Tcoh) + 1;

% Extract SFT coefficients of signal and side bins
S_sig2 = zeros(1, Nseg);
S_side2 = zeros(2*nBinSide, Nseg);
S_sig2_interp = zeros(1, Nseg);
S_side2_interp = zeros(1, Nseg);
for k = 1:Nseg
    S_sig2(k) = S2(fIndex2(k), k);
    S_sig2_interp(k) = S2_interp(fIndex2(k), k);

    for j = 1:nBinSide
        S_side2(nBinSide-j+1, k) = S2(fIndex2(k)-j, k);
        S_side2(nBinSide+j, k) = S2(fIndex2(k)+j, k);
        S_side2_interp(nBinSide-j+1, k) = S2_interp(fIndex2(k)-j, k);
        S_side2_interp(nBinSide+j, k) = S2_interp(fIndex2(k)+j, k);
    end
end
k = 1;
assert(S_sig2_interp(k) == S2_interp(fIndex2(k), k));

% Plot magnitude and phase differences of signal bins
figure;
subplot(1, 3, 1);
p1 = plot(1:Nseg, (abs(S_sig2_interp)./abs(S_sig2)).^2, '-r');
title(sprintf('Signal SFT: %s/%s\n%s = %e%', '$|\tilde{S_2}|^2$', '$|S_2|^2$', '$\dot{f}_2$', fdot2), 'Interpreter', 'Latex');
xlabel('Time Segment', 'Interpreter', 'Latex');
ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
p1.LineWidth = 2;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 18;
grid on;
assert(S_sig2(k) == S2(fIndex2(k), k));

subplot(1, 3, 2);
p2 = plot(1:Nseg, phasediff(S_sig2_interp, S_sig2), '-g');
title(sprintf('Signal SFT: %s - %s\n%s = %e%', '$Phase(\tilde{S_2})$', '$Phase(S_2)$', '$\dot{f}_2$', fdot2), 'Interpreter', 'Latex');
xlabel('Time Segment', 'Interpreter', 'Latex');
ylabel('Phase Difference', 'Interpreter', 'Latex');
p2.LineWidth = 2;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 18;
grid on;

segs = 1:Nseg;
dphi2 = phasediff(S_sig2_interp, S_sig2);

subplot(1, 3, 3);
p3 = plot(segs(abs(dphi2) < phimax), dphi2(abs(dphi2) < phimax), '-b');
title(sprintf('Signal SFT: %s - %s (zoom in)\n%s = %e%\n', '$Phase(\tilde{S_2})$', '$Phase(S_2)$', '$\dot{f}_2$', fdot2), 'Interpreter', 'Latex');
xlabel('Time Segment', 'Interpreter', 'Latex');
ylabel('Phase Difference', 'Interpreter', 'Latex');
p3.LineWidth = 2;
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 18;
grid on;

% Plot magnitude and phase differences of side bins
dphi2_side = cell(1, 2*nBinSide);
figure;
sgtitle(sprintf('Side Bin SFT:\n%s = %e\n blue plots ignore large differences%', '$\dot{f}_2$', fdot2), 'Interpreter', 'Latex');
for j = 1:2*nBinSide
    subplot(3, 2*nBinSide, j);
    p1 = plot(1:Nseg, (abs(S_side2_interp(j, :))./abs(S_side2(j, :))).^2, '-r');
    title(sprintf('Side Bin %d: %s/%s', j, '$|\tilde{S_2}|^2$', '$|S_2|^2$'), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Magnitude Ratio', 'Interpreter', 'Latex');
    p1.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 8;
    grid on;
    assert(S_sig2(k) == S2(fIndex2(k), k));
    
    subplot(3, 2*nBinSide, j+2*nBinSide);
    p2 = plot(1:Nseg, phasediff(S_side2_interp(j, :), S_side2(j, :)), '-g');
    title(sprintf('Side Bin %d: %s - %s', j, '$Phase(\tilde{S_2})$', '$Phase(S_2)$'), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Phase Difference', 'Interpreter', 'Latex');
    p2.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 8;
    grid on;

    dphi2_side{j} = phasediff(S_side2_interp(j, :), S_side2(j, :));
    
    subplot(3, 2*nBinSide, j+4*nBinSide);
    p3 = plot(segs(abs(dphi2_side{j}) < phimax), dphi2_side{j}(abs(dphi2_side{j}) < phimax), '-b');
    title(sprintf('Side Bin %d: %s - %s\n%s = %e%', j, '$Phase(\tilde{S_2})$', '$Phase(S_2)$'), 'Interpreter', 'Latex');
    xlabel('Time Segment', 'Interpreter', 'Latex');
    ylabel('Phase Difference', 'Interpreter', 'Latex');
    p3.LineWidth = 2;
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 8;
    grid on;
end




