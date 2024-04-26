%% Plot and fit maximum time segments that can be interpolated 

clear;
clc;
close all;

Save = 0;

%% Load and Organize Data

N_interp = 5;

Data = readmatrix(sprintf('MaxSegInterp%d.csv', N_interp));
Data(:, 4) = Data(:, 4)-1; % Recorded value was segment of large phase difference. We want the maximum number of safe-to-interpolate segments
fdotvec = unique(Data(:, 1));
Tobsvecs = cell(length(fdotvec), 1);
legends = cell(length(fdotvec), 1);
gvecs = cell(length(fdotvec), 1);
for i = 1:length(fdotvec)
    legends{i} = sprintf('%s = %0.e Hz/s', '$\dot{f}$', fdotvec(i));
    ind1 = (Data(:, 1) == fdotvec(i));
    Tobsvecs{i} = unique(Data(ind1, 2));
    gvecs{i} = cell(length(Tobsvecs{i}), 1);
    for j = 1:length(Tobsvecs{i})
        ind2 = logical((Data(:, 2) == Tobsvecs{i}(j)).*ind1);
        gvecs{i}{j} = Data(ind2, 3);
    end
end

%% Perform LSQ fit and plot

rainbowsc = {0.5*[1, 0, 0], 0.5*[1, 1, 0], 0.5*[0, 1, 0], 0.5*[0, 0, 1], 0.5*[1, 0, 1]};
rainbowsf = {'red', 'yellow', 'green', 'blue', 'magenta'};

LSQ = [log(abs(Data(:, 1:3))), ones(height(Data), 1)];
[param, dparam] = lscov(LSQ, log(Data(:, 4)));

[Tobsmat, gmat] = meshgrid(logspace(log10(1), log10(40), 100), logspace(log10(1), log10(1200), 100));
model_func = @(a, fdot, Tobs_hr, g) (abs(fdot).^a(1)).*(abs(Tobs_hr).^a(2)).*(abs(g).^a(3)).*exp(a(4));

% Plot segments
sc = cell(length(fdotvec), 1);
sf = cell(length(fdotvec), 1);
figure;
hold on;
for i = 1:length(fdotvec)
    ind1 = (Data(:, 1) == fdotvec(i));
    sc{i} = scatter3(Data(ind1, 2), Data(ind1, 3), Data(ind1, 4));
    sf{i} = surf(Tobsmat, gmat, model_func(param, fdotvec(i), Tobsmat, gmat));
    sc{i}.MarkerEdgeColor = rainbowsc{i};
    sc{i}.MarkerFaceColor = rainbowsc{i};
    sf{i}.EdgeAlpha = 0.125;
    sf{i}.FaceAlpha = 0.25;
    sf{i}.FaceColor = rainbowsf{i};
    sc{i}.LineWidth = 3;
    sf{i}.EdgeAlpha = 0.5;
    sf{i}.EdgeColor = rainbowsf{i};
end
hold off;
legend([sf{1}, sf{2}, sf{3}, sf{4}, sf{5}], legends, 'Interpreter', 'LaTex', 'Location', 'eastoutside');
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.ZScale = 'log';
ax.LineWidth = 3;
ax.FontSize = 18;
title(sprintf('Maximum number of segments that can be interpolated\n without phase inaccuracies'), 'Interpreter', 'Latex');
xlabel('$T_{obs}$ (hr)', 'Interpreter', 'Latex');
ylabel('$g$', 'Interpreter', 'Latex');
zlabel('Segments', 'Interpreter', 'Latex');
grid on;

% Plot time
sc = cell(length(fdotvec), 1);
sf = cell(length(fdotvec), 1);
figure;
hold on;
for i = 1:length(fdotvec)
    ind1 = (Data(:, 1) == fdotvec(i));
    sc{i} = scatter3(Data(ind1, 2), Data(ind1, 3), Data(ind1, 4)/sqrt(2*abs(fdotvec(i)))/3600);
    sf{i} = surf(Tobsmat, gmat, model_func(param, fdotvec(i), Tobsmat, gmat)/sqrt(2*abs(fdotvec(i)))/3600);
    sc{i}.MarkerEdgeColor = rainbowsc{i};
    sc{i}.MarkerFaceColor = rainbowsc{i};
    sf{i}.EdgeAlpha = 0.125;
    sf{i}.FaceAlpha = 0.25;
    sf{i}.FaceColor = rainbowsf{i};
    sc{i}.LineWidth = 3;
    sf{i}.EdgeAlpha = 0.5;
    sf{i}.EdgeColor = rainbowsf{i};
end
sf2 = surf(Tobsmat, gmat, Tobsmat);
sf2.EdgeAlpha = 0.25;
sf2.FaceAlpha = 0.25;
sf2.FaceColor = 'black';
hold off;
legend([sf{1}, sf{2}, sf{3}, sf{4}, sf{5}, sf2], {legends{1:end}, '$T_{obs, max} = T_{obs}$'}, 'Interpreter', 'LaTex', 'Location', 'eastoutside');
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.ZScale = 'log';
ax.LineWidth = 3;
ax.FontSize = 18;
title(sprintf('Maximum observation time that can be interpolated\n without phase inaccuracies'), 'Interpreter', 'Latex');
xlabel('$T_{obs}$ (hr)', 'Interpreter', 'Latex');
ylabel('$g$', 'Interpreter', 'Latex');
zlabel('$T_{obs, max}$ (hr)', 'Interpreter', 'Latex');
grid on;

if (Save)
    fout = fopen('MaxSegInterpParams.csv', 'a');
    fprintf(fout, '%d, %e, %e, %e, %e\n', N_interp, param(1), param(2), param(3), param(4));
end
