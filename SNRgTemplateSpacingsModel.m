% Takes Template Spacings and SNRgs and fits them to model

clear;
clc;
close all;

ParsevalSNR = 1;
PlotChiSqrContributions = 0;
CrossTerms = 0;

%% Fit Template Spacings

% Read in and configure data
fdotvec_sig = -5.*10.^(-9:-5);
Ntrials = 10;
X = [];
Y = [];
dY = [];
Tobsvecs = cell(length(fdotvec_sig));
Nsegvecs = Tobsvecs;
for j = 1:length(fdotvec_sig)
    fdot_sig = fdotvec_sig(j);
    Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
    Tcoh = Tcoh_hr * 3600.;
    Tcoh = floor(Tcoh);
    Tcoh_hr = Tcoh/3600.;
    if (ParsevalSNR == 0)
        filenameIn = sprintf('SNRgTemplateSpacingfdot_%0.e.csv', fdot_sig);
    else
        filenameIn = sprintf('SNRgTemplateSpacingfdot_%0.eP.csv', fdot_sig);
    end

    Data = readmatrix(filenameIn);
    Tobsvecs{j} = unique(Data(:, 1));
    Nsegvecs{j} = round(Tobsvecs{j}./Tcoh_hr);
    for k = 1:length(Tobsvecs{j})
        for i = 1:Nsegvecs{j}(k)
            X = [X; log10(abs(fdot_sig)), log10(Tobsvecs{j}(k)), log10(i)];
            dataIndices = logical((Data(:, 2)==i).*(Data(:, 1)==Tobsvecs{j}(k)));
            assert(sum(dataIndices)==Ntrials, ['Error: Incorrect number of trials for Tobs = ', num2str(Tobsvecs{j}(k)), ' hr and i = ', num2str(i)]);
            Y = [Y; log10(mean(Data(dataIndices, 3)))];
            dY = [dY; std(Data(dataIndices, 3))/mean(Data(dataIndices, 3))/log(10)/sqrt(Ntrials)];
        end
    end
end

% Do additional fits to add to uncertainties
plane_model = @(a, X) a(1).*X(:, 1) + a(2).*X(:, 2) + a(3);
plane_fitParams = zeros(length(fdotvec_sig), 3);
for j = 1:length(fdotvec_sig)
    dataIndices = logical(X(:, 1) == log10(abs(fdotvec_sig(j))));
    plane_fitParams(j, :) = nlinfit(X(dataIndices, 2:3), Y(dataIndices), plane_model, -2*rand(1, 3));
    for k = 1:length(Tobsvecs{j})
        for i = 1:Nsegvecs{j}(k)
            dataIndex = logical(dataIndices .* (X(:, 3) == log10(i)) .* (X(:, 2) == log10(Tobsvecs{j}(k))));
            sysSigma = abs(Y(dataIndex) - plane_model(plane_fitParams(j, :), log10([Tobsvecs{j}(k), i])));
            dY(dataIndex) = sqrt(dY(dataIndex)^2 + sysSigma^2);
        end
    end
end

% Define model function and perform chi-squared minimization
lin_model = @(a, X) a(1) + a(2).*X(:);
initialFitParams = zeros(1, 4+2*CrossTerms);
for l = 1:3
    initialFitParams((2*l-1):(2*l)) = nlinfit(log10(abs(fdotvec_sig)), plane_fitParams(:, l), lin_model, -rand(1, 2));
end
if(CrossTerms == 0)
    initialFitParams = initialFitParams([1 3 5 6]);
end

model_func = @(a, X) a(1).*X(:, 2) + ...
                     a(2).*X(:, 3) + ...
                     a(3) + a(4).*X(:, 1);

[fitParams_TS, dParams_TS] = fitChiSquare(X, Y, model_func, initialFitParams, zeros(height(X), 3), dY);


% Plot data points and fits on surfaces of constant fdot
if (CrossTerms == 1)
    model_func2 = @(a, fdot_log, Tobs_log, i_log) (a(1) + a(2).*fdot_log).*Tobs_log + ...
                                                  (a(3) + a(4).*fdot_log).*i_log + ...
                                                  (a(5) + a(6).*fdot_log);
else
    model_func2 = @(a, fdot_log, Tobs_log, i_log) a(1).*Tobs_log + ...
                                                  a(2).*i_log + ...
                                                  a(3) + a(4).*fdot_log;
end

modelmats = cell(length(fdotvec_sig));
for k = 1:length(fdotvec_sig)
    [Tobsmat, imat] = meshgrid(linspace(log10(1), log10(max(Tobsvecs{k})), 100), linspace(log10(1), log10(max(Nsegvecs{k})), 100));
    modelmats{k} = model_func2(fitParams_TS, log10(abs(fdotvec_sig(k))), Tobsmat, imat);
    
    figure;
    s1000 = scatter3(X(X(:, 1)==log10(abs(fdotvec_sig(k))), 2), X(X(:, 1)==log10(abs(fdotvec_sig(k))), 3), Y(X(:, 1)==log10(abs(fdotvec_sig(k)))), '*m');
    hold on;
    sf1000 = surf(Tobsmat, imat, modelmats{k});
    title(['Template Spacing over Tobs and i for fdot = -5.e', num2str(log10(abs(fdotvec_sig(k))/5)), ' Hz/s']);
    xlabel('log(Tobs) (log(hr))');
    ylabel('log(i)');
    zlabel('log(dfdot) (log(Hz/s))');
    grid on;
    s1000.LineWidth = 3;
    sf1000.EdgeAlpha = 0.5;
    ax = gca; 
    ax.LineWidth = 3;
    ax.FontSize = 16;
    hold off;
end

% Calculate chi^2 contributions
if (PlotChiSqrContributions == 1)
    chisqrContributions = cell(1, length(fdotvec_sig));
    for j = 1:length(fdotvec_sig)
        chisqrContributions{j} = zeros(length(Tobsvecs{j}), Nsegvecs{j}(end));
        for k = 1:length(Tobsvecs{j})
            for i = 1:Nsegvecs{j}(k)
                dataIndex = logical( (X(:, 1) == log10(abs(fdotvec_sig(j)))) .* ...
                                     (X(:, 2) == log10(Tobsvecs{j}(k))) .* ...
                                     (X(:, 3) == log10(i)) );
                chisqrContributions{j}(k, i) = ((Y(dataIndex) - model_func(fitParams_TS, X(dataIndex, :)))/dY(dataIndex))^2;
            end
        end
    end
    
    for j = 1:length(fdotvec_sig)
        [imat, Tobsmat] = meshgrid(1:Nsegvecs{j}(end), Tobsvecs{j});
        figure;
        sfk = surf(Tobsmat, imat, chisqrContributions{j});
        title(sprintf('Chi^2 Contributions from Each Data Point fdot = %2.e Hz/s', fdotvec_sig(j)));
        xlabel('Tobs (hr)');
        ylabel('i')
        ax = gca;
        ax.FontSize = 16;
        ax.LineWidth = 3;
        grid on;
    end
end

%% Fit SNRg Values

% Read in and configure data
X = [];
Y = [];
dY = [];
for j = 1:length(fdotvec_sig)
    fdot_sig = fdotvec_sig(j);
    Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
    Tcoh = Tcoh_hr * 3600.;
    Tcoh = floor(Tcoh);
    Tcoh_hr = Tcoh/3600.;
    for k = 1:length(Tobsvecs{j})
        if (ParsevalSNR == 1)
            filenameIn = sprintf('SNRgfdot_%0.es-2Tobs_%0.fhrP.csv', fdotvec_sig(j), Tobsvecs{j}(k));
        else
            filenameIn = sprintf('SNRgfdot_%0.es-2Tobs_%0.fhr.csv', fdotvec_sig(j), Tobsvecs{j}(k));
        end
        Data = readmatrix(filenameIn);
        assert(height(Data)==Ntrials, ['Error: Incorrect number of trials for Tobs = ', num2str(Tobsvecs{j}(k)), ' hr\n']);
        for i = 1:Nsegvecs{j}(k)
            X = [X; log10(abs(fdot_sig)), log10(Tobsvecs{j}(k)), log10(i)];
            Y = [Y; log10(mean(Data(:, i)))];
            dY = [dY; std(Data(:, i))/mean(Data(:, i))/log(10)/sqrt(Ntrials)];
        end
    end
end

% Do additional fits to add to uncertainties
plane_model = @(a, X) a(1).*X(:, 1) + a(2).*X(:, 2) + a(3);
plane_fitParams = zeros(length(fdotvec_sig), 3);
for j = 1:length(fdotvec_sig)
    dataIndices = logical(X(:, 1) == log10(abs(fdotvec_sig(j))));
    plane_fitParams(j, :) = nlinfit(X(dataIndices, 2:3), Y(dataIndices), plane_model, -2*rand(1, 3));
    for k = 1:length(Tobsvecs{j})
        for i = 1:Nsegvecs{j}(k)
            dataIndex = logical(dataIndices .* (X(:, 3) == log10(i)) .* (X(:, 2) == log10(Tobsvecs{j}(k))));
            sysSigma = abs(Y(dataIndex) - plane_model(plane_fitParams(j, :), log10([Tobsvecs{j}(k), i])));
            dY(dataIndex) = sqrt(dY(dataIndex)^2 + sysSigma^2);
        end
    end
end

% Define model function and perform chi-squared minimization
lin_model = @(a, X) a(1) + a(2).*X(:);
initialFitParams = zeros(1, 4+2*CrossTerms);
for l = 1:3
    initialFitParams((2*l-1):(2*l)) = nlinfit(log10(abs(fdotvec_sig)), plane_fitParams(:, l), lin_model, -rand(1, 2));
end
if(CrossTerms == 0)
    initialFitParams = initialFitParams([1 3 5 6]);
end

model_func = @(a, X) a(1).*X(:, 2) + ...
                     a(2).*X(:, 3) + ...
                     a(3) + a(4).*X(:, 1);

[fitParams_SNR, dParams_SNR] = fitChiSquare(X, Y, model_func, initialFitParams, zeros(height(X), 3), dY);


% Plot data points and fits on surfaces of constant fdot
if (CrossTerms == 1)
    model_func2 = @(a, fdot_log, Tobs_log, i_log) (a(1) + a(2).*fdot_log).*Tobs_log + ...
                                                  (a(3) + a(4).*fdot_log).*i_log + ...
                                                  (a(5) + a(6).*fdot_log);
else
    model_func2 = @(a, fdot_log, Tobs_log, i_log) a(1).*Tobs_log + ...
                                                  a(2).*i_log + ...
                                                  a(3) + a(4).*fdot_log;
end

modelmats = cell(length(fdotvec_sig));
for k = 1:length(fdotvec_sig)
    [Tobsmat, imat] = meshgrid(linspace(log10(1), log10(max(Tobsvecs{k})), 100), linspace(log10(1), log10(max(Nsegvecs{k})), 100));
    modelmats{k} = model_func2(fitParams_SNR, log10(abs(fdotvec_sig(k))), Tobsmat, imat);
    
    figure;
    s1000 = scatter3(X(X(:, 1)==log10(abs(fdotvec_sig(k))), 2), X(X(:, 1)==log10(abs(fdotvec_sig(k))), 3), Y(X(:, 1)==log10(abs(fdotvec_sig(k)))), '*m');
    hold on;
    sf1000 = surf(Tobsmat, imat, modelmats{k});
    title(['SNRg_i over Tobs and i for fdot = -5.e', num2str(log10(abs(fdotvec_sig(k))/5)), ' Hz/s']);
    xlabel('log(Tobs) (log(hr))');
    ylabel('log(i)');
    zlabel('log(dfdot) (log(Hz/s))');
    grid on;
    s1000.LineWidth = 3;
    sf1000.EdgeAlpha = 0.5;
    ax = gca; 
    ax.LineWidth = 3;
    ax.FontSize = 16;
    hold off;
end

% Calculate chi^2 contributions
if (PlotChiSqrContributions == 1)
    chisqrContributions = cell(1, length(fdotvec_sig));
    for j = 1:length(fdotvec_sig)
        chisqrContributions{j} = zeros(length(Tobsvecs{j}), Nsegvecs{j}(end));
        for k = 1:length(Tobsvecs{j})
            for i = 1:Nsegvecs{j}(k)
                dataIndex = logical( (X(:, 1) == log10(abs(fdotvec_sig(j)))) .* ...
                                     (X(:, 2) == log10(Tobsvecs{j}(k))) .* ...
                                     (X(:, 3) == log10(i)) );
                chisqrContributions{j}(k, i) = ((Y(dataIndex) - model_func(fitParams_SNR, X(dataIndex, :)))/dY(dataIndex))^2;
            end
        end
    end
    
    for j = 1:length(fdotvec_sig)
        [imat, Tobsmat] = meshgrid(1:Nsegvecs{j}(end), Tobsvecs{j});
        figure;
        sfk = surf(Tobsmat, imat, chisqrContributions{j});
        title(sprintf('Chi^2 Contributions from Each Data Point fdot = %2.e Hz/s', fdotvec_sig(j)));
        xlabel('Tobs (hr)');
        ylabel('i')
        ax = gca;
        ax.FontSize = 16;
        ax.LineWidth = 3;
        grid on;
    end
end

