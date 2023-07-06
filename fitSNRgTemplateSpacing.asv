%% SNRg Template Spacing Fits

clear;
clc;
close all;

ParsevalSNR = 1;
Fits2D = 0;
RunningAveSigma = 0;
IterFitSigma = 1;
OutputFile = 0;

% Choose exactly one of the following fit algorithms
UsefitChiSquare = 1;
Usenlinfit = 0;
UseLSQCoeffs = 0;

% Choose degree of polynomial surface fit (1, 2, 3, 
% pi for 2 piecewise planes intersecting at constant i
% 2*pi for 2 piecewise planes intersecting at arbitrary line parallel to i-Tobs plane)
fitDegree = 2*pi;

%% Read in and configure data

% Select fdot value (Hz/s)
fdot_sig = -5.e-8; % One of the following: -5.e-9, -5.e-8, -5.e-7, -5.e-6, -5.e-5

Tcoh_hr = 1./3600.*sqrt(0.5/abs(fdot_sig));
Tcoh = Tcoh_hr * 3600.;
Tcoh = floor(Tcoh);
Tcoh_hr = Tcoh/3600.;

% Select number of trials used and arrange data
Ntrials = 10;
if (ParsevalSNR == 0)
    filenameIn = sprintf('SNRgTemplateSpacingfdot_%0.e.csv', fdot_sig);
else
    filenameIn = sprintf('SNRgTemplateSpacingfdot_%0.eP.csv', fdot_sig);
end
Data = readmatrix(filenameIn);
Tobsvec = unique(Data(:, 1));
Nsegvec = round(Tobsvec./Tcoh_hr);

% Calculate uncertainties and means
uncertainties = zeros(length(Tobsvec), max(Nsegvec));
data = zeros(length(Tobsvec), max(Nsegvec));
for k = 1:length(Tobsvec)
    for i = 1:Nsegvec(k)
        dataIndices = logical((Data(:, 2)==i).*(Data(:, 1)==Tobsvec(k)));
        assert(sum(dataIndices)==Ntrials, ['Error: Incorrect number of trials for Tobs = ', num2str(Tobsvec(k)), ' hr and i = ', num2str(i)]);
        uncertainties(k, i) = std(Data(dataIndices, 3));
        data(k, i) = mean(Data(dataIndices, 3));
    end
end
uncertainties = uncertainties./sqrt(Ntrials);

%% Calculate Other Uncertainties

if (RunningAveSigma == 1)
    runningAveSigma = zeros(size(uncertainties));
    for k = 1:length(Tobsvec)
        dfdotave = zeros(1, round(Nsegvec(k)));
        for i = 1:Nsegvec(k)
            if (i == 1 || i == Nsegvec(k))
                dfdotave(i) = log10(data(k, i));
            elseif (i == 2 || i == (Nsegvec(k)-1))
                dfdotave(i) = mean(log10(data(k, (i-1):(i+1))));
            elseif (i == 3 || i == (Nsegvec(k)-2))
                dfdotave(i) = mean(log10(data(k, (i-2):(i+2))));
            elseif (i == 4 || i == (Nsegvec(k)-3))
                dfdotave(i) = mean(log10(data(k, (i-3):(i+3))));
            else
                dfdotave(i) = mean(log10(data(k, (i-4):(i+4))));
            end
        end
        
        figure(1000 + k)
        p1000k = plot(log10(1:Nsegvec(k)), dfdotave, '--r');
        hold on; 
        p1000k2 = plot(log10(1:Nsegvec(k)), log10(data(k, 1:Nsegvec(k))), '-k');
        p1000k2.LineWidth = 2;
        title(['Running Average of Data for Tobs = ', num2str(Tobsvec(k)), ' hr']);
        xlabel('Base 10 log of i');
        ylabel('Base 10 log of Template Spacing (log(Hz/s))');
        legend('Running Average', 'Data');
        grid on;
        hold off;
        
        runningAveSigma(k, 1:Nsegvec(k)) = abs(log10(data(k, 1:Nsegvec(k)))-dfdotave);
        runningAveSigma(k, 1) = runningAveSigma(k, 2);
        runningAveSigma(k, Nsegvec(k)) = runningAveSigma(k, Nsegvec(k)-1);
        
    end
end

if (IterFitSigma == 1)
    iterFitSigma = zeros(size(uncertainties));
    quad_model = @(a, x) a(1)*x.^2 + a(2)*x + a(3);
    for k = 1:length(Tobsvec)
        errorfitparams = nlinfit(log10(1:Nsegvec(k)), log10(data(k, 1:Nsegvec(k))), quad_model, 2.*rand(1, 3)-1);        
        
        figure(1000 + k)
        p1000k = plot(log10(1:Nsegvec(k)), quad_model(errorfitparams, log10(1:Nsegvec(k))), '--r');
        hold on; 
        p1000k2 = plot(log10(1:Nsegvec(k)), log10(data(k, 1:Nsegvec(k))), '-k');
        p1000k2.LineWidth = 2;
        title(['Quadratic Fit of Template Spacings for Tobs = ', num2str(Tobsvec(k)), ' hr']);
        xlabel('Base 10 log of i');
        ylabel('Base 10 log of Template Spacing (log(Hz/s))');
        legend('Quadratic fit', 'Data');
        grid on;
        hold off;
        
        iterFitSigma(k, 1:Nsegvec(k)) = abs(log10(data(k, 1:Nsegvec(k)))-quad_model(errorfitparams, log10(1:Nsegvec(k))));
    end
end

%% Perform 2D fits

if (Fits2D == 1)
    Nfits2D = 9;
    fitCoeffs2Di = cell(1, round(Nsegvec(1)));
    fitCoeffs2DTobs = cell(1, length(Tobsvec));
    model_func2 = @(a, x) a(1).*x + a(2);
    for i = 1:Nsegvec(1)
        dfdot = data(:, i);
        uncertainty = uncertainties(:, i);
        Tobs = Tobsvec(dfdot ~= 0);
        uncertainty = uncertainty(dfdot ~= 0);
        dfdot = dfdot(dfdot ~= 0);

        [fitCoeffs, dParams] = fitChiSquare(log10(Tobs), log10(dfdot), model_func2, 2.*rand(1, 2)-1, zeros(size(Tobs)), uncertainty./dfdot./log(10));
        dslope = dParams(1).d;
        doffset = dParams(2).d;

        fitCoeffs2Di{i} = fitCoeffs;

        figure(i)
        ei = errorbar(log10(Tobs), log10(dfdot), uncertainty./dfdot./log(10));
        hold on;
        pi = plot(log10(Tobs), fitCoeffs(1).*log10(Tobs) + fitCoeffs(2), '--r');
        title(['Template Spacing vs Tobs for i = ', num2str(i)]);
        xlabel('log(Tobs) (log(hr))');
        ylabel('log(dfdot) (log(Hz/s))');
        grid on;
        si.LineWidth = 3;
        pi.LineWidth = 1;
        legend('Calculated Data', ['Fit (slope = ', num2str(fitCoeffs(1)), '+/-', num2str(dslope), ')']);
        ax = gca; 
        ax.LineWidth = 3;
        ax.FontSize = 16;
    end
    
    for k = 1:length(Tobsvec)
        dfdot = data(k, :);
        dfdot = dfdot(dfdot ~= 0);
        uncertainty = uncertainties(k, 1:length(dfdot));
        ivec = 1:length(dfdot);
        
        [fitCoeffs, dParams] = fitChiSquare(log10(ivec), log10(dfdot), model_func2, 2.*rand(1, 2)-1, zeros(size(ivec)), uncertainty./dfdot./log(10));
        dslope = dParams(1).d;
        doffset = dParams(2).d;

        fitCoeffs2DTobs{i} = fitCoeffs;

        figure(k + round(Nsegvec(1)))
        ek = errorbar(log10(ivec), log10(dfdot), uncertainty./dfdot./log(10));
        hold on;
        pk = plot(log10(ivec), fitCoeffs(1).*log10(ivec) + fitCoeffs(2), '--r');
        title(['Template Spacing vs i for Tobs = ', num2str(Tobsvec(k))]);
        xlabel('log(i)');
        ylabel('log(dfdot) (log(Hz/s))');
        grid on;
        pk.LineWidth = 1;
        legend('Calculated Data', ['Fit (slope = ', num2str(fitCoeffs(1)), '+/-', num2str(dslope), ')']);
        ax = gca; 
        ax.LineWidth = 3;
        ax.FontSize = 16;
    end
end


%% Perform Full Fit

% Reformat data into vectors for fitting
X = [];
Y = [];
dY = [];
for m = 1:length(data(:, 1))
    for n = 2:length(data(1, :))
        if (n <= round(Nsegvec(m)))
            X = [X; log10(Tobsvec(m)), log10(n)];
            Y = [Y; log10(data(m, n))];
            if (RunningAveSigma == 1)
                dY = [dY; runningAveSigma(m, n)];
            elseif (IterFitSigma == 1)
                dY = [dY; sqrt((uncertainties(m, n)./data(m, n)/log(10))^2 + iterFitSigma(m, n)^2)];
            else
                dY = [dY; uncertainties(m, n)./data(m, n)/log(10)];
            end
        end
    end
end

% Set model function and least squares fit based on degree of fit
if (fitDegree == 3)
    lsqCoeffs = [X(:, 1).^3, (X(:, 1).^2).*X(:, 2), X(:, 1).*(X(:, 2).^2), X(:, 2).^3, X(:, 1).^2, X(:, 1).*X(:, 2), X(:, 2).^2, X(:, 1), X(:, 2), ones(size(X(:, 2)))] \ Y;

    model_func = @(a, X) a(1).*X(:, 1).^3 + a(2).*(X(:, 1).^2).*X(:, 2) + a(3).*X(:, 1).*(X(:, 2).^2) + a(4).*X(:, 2).^3 + a(5).*X(:, 1).^2 + a(6).*X(:, 1).*X(:, 2) + a(7).*X(:, 2).^2 + a(8).*X(:, 1) + a(9).*X(:, 2) + a(10);
    d = zeros(1, 10);
elseif (fitDegree == 2)
    lsqCoeffs = [X(:, 1).^2, X(:, 1).*X(:, 2), X(:, 2).^2, X(:, 1), X(:, 2), ones(size(X(:, 2)))] \ Y;

    model_func = @(a, X)  a(1).*X(:, 1).^2 + a(2).*X(:, 1).*X(:, 2) + a(3).*X(:, 2).^2 + a(4).*X(:, 1) + a(5).*X(:, 2) + a(6);
    d = zeros(1, 6);
elseif (fitDegree == 1)
    lsqCoeffs = [X(:, 1), X(:, 2), ones(size(X(:, 2)))] \ Y;

    model_func = @(a, X)  a(1).*X(:, 1) + a(2).*X(:, 2) + a(3);
    d = zeros(1, 3);
elseif (fitDegree == pi)
    lsqCoeffs = [X(:, 1), X(:, 2), ones(size(X(:, 2)))] \ Y;
    
    if (fdot_sig == -5.e-5)
        lsqCoeffs = [lsqCoeffs; 2.; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-6)
        lsqCoeffs = [lsqCoeffs; 1.5; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-7)
        lsqCoeffs = [lsqCoeffs; 1.; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-8)
        lsqCoeffs = [lsqCoeffs; 0.75; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-9)
        lsqCoeffs = [lsqCoeffs; 0.5; lsqCoeffs(2)];
    end
    
    model_func = @(a, X) (a(1).*X(:, 1) + a(2).*X(:, 2) + a(3)).*(X(:, 2) < a(4)) + (a(1).*X(:, 1) + a(5).*X(:, 2) + a(2)*a(4) + a(3) - a(5)*a(4)).*(X(:, 2) >= a(4));
    d = zeros(1, 5);
elseif (fitDegree == 2*pi)
    lsqCoeffs = [X(:, 1), X(:, 2), ones(size(X(:, 2)))] \ Y;

    if (fdot_sig == -5.e-5)
        lsqCoeffs = [lsqCoeffs; 0; 2.090003; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-6)
        lsqCoeffs = [lsqCoeffs; 0; 1.572379; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-7)
        lsqCoeffs = [lsqCoeffs; 0; 1.097328; lsqCoeffs(2)];
    elseif (fdot_sig == -5.e-8)
        lsqCoeffs = [lsqCoeffs; 0; 0.674366; lsqCoeffs(2)];
    end

    model_func = @(a, X) (a(1).*X(:, 1) + a(2).*X(:, 2) + a(3)) .* ( X(:, 2) < a(4).*X(:, 1) + a(5) ) + ((a(1) + a(2)*a(4) - a(6)*a(4)).*X(:, 1) + a(6).*X(:, 2) + a(2)*a(5) + a(3) - a(6)*a(5)) .* ( X(:, 2) >= a(4).*X(:, 1) + a(5) );
    d = zeros(1, 6);
elseif (fitDegree == 3*pi)
    lsqCoeffs = [X(:, 2), ones(size(X(:, 2)))] \ Y;

    if (fdot_sig == -5.e-5)
        lsqCoeffs = [lsqCoeffs; 2.; -1.5];
    elseif (fdot_sig == -5.e-6)
        lsqCoeffs = [lsqCoeffs; 1.5; -2];
    elseif (fdot_sig == -5.e-7)
        lsqCoeffs = [lsqCoeffs; 1.; -1.5];
    elseif (fdot_sig == -5.e-8)
        lsqCoeffs = [lsqCoeffs; 0.75; -1.5];
    elseif (fdot_sig == -5.e-9)
        lsqCoeffs = [lsqCoeffs; 0.8; -1.5];
    end

    model_func = @(a, X) (a(1).*X(:, 2) + a(2)).*(X(:, 2) < a(3)) + (a(4).*X(:, 2) + a(1)*a(3) + a(2) - a(4)*a(3)).*(X(:, 2) >= a(3));
    d = zeros(1, 4);
end

% Perform fit
if (Usenlinfit == 1)
    [fitCoeffs, R, J, Cov, MSE, ErrorModelInfo] = nlinfit(X, Y, model_func, lsqCoeffs, 'Weights', dY^(-1));
    for k = 1:length(d)
        d(k) = sqrt(Cov(k,k));
    end
elseif (UsefitChiSquare == 1)
    [fitCoeffs, dParams] = fitChiSquare(X, Y, model_func, lsqCoeffs, zeros(size(X)), dY);
    for k = 1:length(d)
        d(k) = dParams(k).d;
    end
elseif (UseLSQCoeffs == 1)
    fitCoeffs = lsqCoeffs;
end

% Calculate surface based on fit parameters
[Tobsmat, imat] = meshgrid(linspace(log10(min(Tobsvec)), log10(max(Tobsvec)), 100), linspace(log10(1), log10(max(Nsegvec)), 100));
if (fitDegree == 3)
    model_func2 = @(a, X, Y) a(1).*X.^3 + a(2).*(X.^2).*Y + a(3).*X.*(Y.^2) + a(4).*Y.^3 + a(5).*X.^2 + a(6).*X.*Y + a(7).*Y.^2 + a(8).*X + a(9).*Y + a(10);
elseif (fitDegree == 2)
    model_func2 = @(a, X, Y)  a(1).*X.^2 + a(2).*X.*Y + a(3).*Y.^2 + a(4).*X + a(5).*Y + a(6);
elseif (fitDegree == 1)
    model_func2 = @(a, X, Y)  a(1).*X + a(2).*Y + a(3);
elseif (fitDegree == pi)
    model_func2 = @(a, X, Y) (a(1).*X + a(2).*Y + a(3)).*(Y < a(4)) + (a(1).*X + a(5).*Y + a(2)*a(4)+a(3)-a(5)*a(4)).*(Y >= a(4));
elseif (fitDegree == 2*pi)
    model_func2 = @(a, X, Y) (a(1).*X + a(2).*Y + a(3)) .* ( Y < a(4).*X + a(5) ) + ((a(1) + a(2)*a(4) - a(6)*a(4)).*X + a(6).*Y + a(2)*a(5) + a(3) - a(6)*a(5)) .* ( Y >= a(4).*X + a(5) );
elseif (fitDegree == 3*pi)
    model_func2 = @(a, X, Y) (a(1).*Y + a(2)).*(Y < a(3)) + (a(4).*Y + a(1)*a(3)+a(2)-a(4)*a(3)).*(Y >= a(3));
end
modelmat = model_func2(fitCoeffs, Tobsmat, imat);

%% Plot Surface 

figure(1000)
s1000 = scatter3(X(:, 1), X(:, 2), Y, '*m');
hold on;
sf1000 = surf(Tobsmat, imat, modelmat);
title('Template Spacing over Tobs and i');
xlabel('log(Tobs) (log(hr))');
ylabel('log(i)');
zlabel('log(dfdot) (log(Hz/s))');
if (fitDegree == pi)
    Tpoints = linspace(log10(min(Tobsvec)), log10(max(Tobsvec)), 100);
    p1000 = plot3(Tpoints, fitCoeffs(4).*ones(1, 100), model_func2(fitCoeffs, Tpoints, fitCoeffs(4).*ones(1, 100)), '-y');
    p1000.LineWidth = 2;
elseif (fitDegree == 2*pi)
    Tpoints = linspace(log10(min(Tobsvec)), log10(max(Tobsvec)), 100);
    p1000 = plot3(Tpoints, fitCoeffs(4).*Tpoints + fitCoeffs(5), model_func2(fitCoeffs, Tpoints, fitCoeffs(4).*Tpoints + fitCoeffs(5)), '-y');
    p1000.LineWidth = 2;
elseif (fitDegree == 3*pi)
    Tpoints = linspace(log10(min(Tobsvec)), log10(max(Tobsvec)), 100);
    p1000 = plot3(Tpoints, fitCoeffs(3).*ones(1, 100), model_func2(fitCoeffs, Tpoints, fitCoeffs(3).*ones(1, 100)), '-y');
    p1000.LineWidth = 2;
end
grid on;
s1000.LineWidth = 3;
sf1000.EdgeAlpha = 0.5;

ax = gca; 
ax.LineWidth = 3;
ax.FontSize = 16;

if (OutputFile == 1)
    if (ParsevalSNR == 1)
        filenameOut = sprintf('TemplateSpacingFitCoeffsfitDegree_%fP.csv', fitDegree);
    else 
        filenameOut = sprintf('TemplateSpacingFitCoeffsfitDegree_%f.csv', fitDegree);
    end
    fout = fopen(filenameOut, 'a');
    fprintf(fout, '%e, ', fdot_sig);
end

% Print model function with optimized parameters
fprintf('Surface fit: log(dfdot) = \n');
if (fitDegree == 3)
    fprintf('\t  (%f +/- %f)*log(Tobs)^3 +\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(Tobs)^2*log(i)\n', fitCoeffs(2), d(2));
    fprintf('\t+ (%f +/- %f)*log(Tobs)*log(i)^2\n', fitCoeffs(3), d(3));
    fprintf('\t+ (%f +/- %f)*log(i)^3\n', fitCoeffs(4), d(4));
    fprintf('\t+ (%f +/- %f)*log(Tobs)^2\n', fitCoeffs(5), d(5));
    fprintf('\t+ (%f +/- %f)*log(Tobs)*log(i)\n', fitCoeffs(6), d(6));
    fprintf('\t+ (%f +/- %f)*log(i)^2\n', fitCoeffs(7), d(7));
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(8), d(8));
    fprintf('\t+ (%f +/- %f)*log(i)\n', fitCoeffs(9), d(9));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(10), d(10));
elseif (fitDegree == 2)
    fprintf('\t+ (%f +/- %f)*log(Tobs)^2\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(Tobs)*log(i)\n', fitCoeffs(2), d(2));
    fprintf('\t+ (%f +/- %f)*log(i)^2\n', fitCoeffs(3), d(3));
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(4), d(4));
    fprintf('\t+ (%f +/- %f)*log(i)\n', fitCoeffs(5), d(5));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(6), d(6));
elseif (fitDegree == 1)
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(i)\n', fitCoeffs(2), d(2));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(3), d(3));
elseif (fitDegree == pi)
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(i)\t}log(i) < (%f +/- %f)\n', fitCoeffs(2), d(2), fitCoeffs(4), d(4));
    fprintf('\t+ (%f +/- %f)\n\n', fitCoeffs(3), d(3));
    
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(i)\t}log(i) >= (%f +/- %f)\n', fitCoeffs(5), d(5), fitCoeffs(4), d(4));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(2)*fitCoeffs(4) + fitCoeffs(3) - fitCoeffs(5)*fitCoeffs(4), sqrt((fitCoeffs(2)*d(4))^2 + (fitCoeffs(4)*d(2))^2 + d(3)^2 + (fitCoeffs(5)*d(4))^2 + (fitCoeffs(4)*d(5))^2));
elseif (fitDegree == 2*pi)
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(1), d(1));
    fprintf('\t+ (%f +/- %f)*log(i)\t}log(i) < (%f +/- %f)*log(Tobs) + (%f +/- %f)\n', fitCoeffs(2), d(2), fitCoeffs(4), d(4), fitCoeffs(5), d(5));
    fprintf('\t+ (%f +/- %f)\n\n', fitCoeffs(3), d(3));
    
    fprintf('\t+ (%f +/- %f)*log(Tobs)\n', fitCoeffs(1) + fitCoeffs(2)*fitCoeffs(4) - fitCoeffs(6)*fitCoeffs(4), sqrt(d(1)^2 + (fitCoeffs(2)*d(4))^2 + (fitCoeffs(4)*d(2))^2 + (fitCoeffs(6)*d(4))^2 + (fitCoeffs(4)*d(6))^2));
    fprintf('\t+ (%f +/- %f)*log(i)\t}log(i) >= (%f +/- %f)*log(Tobs) + (%f +/- %f)\n', fitCoeffs(6), d(6),  fitCoeffs(4), d(4), fitCoeffs(5), d(5));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(2)*fitCoeffs(5) + fitCoeffs(3) - fitCoeffs(6)*fitCoeffs(5), sqrt((fitCoeffs(2)*d(5))^2 + (fitCoeffs(5)*d(2))^2 + d(3)^2 + (fitCoeffs(6)*d(5))^2 + (fitCoeffs(5)*d(6))^2));
elseif (fitDegree == 3*pi)
    fprintf('\t(%f +/- %f)*log(i)\n', fitCoeffs(1), d(1));
    fprintf('\t\t\t\t\t } log(i) < (%f +/- %f)\n', fitCoeffs(3), d(3));
    fprintf('\t+ (%f +/- %f)\n\n', fitCoeffs(2), d(2));
    
    fprintf('\t(%f +/- %f)*log(i)\n', fitCoeffs(4), d(4));
    fprintf('\t\t\t\t\t } log(i) >= (%f +/- %f)\n', fitCoeffs(3), d(3));
    fprintf('\t+ (%f +/- %f)\n', fitCoeffs(1)*fitCoeffs(3) + fitCoeffs(2) - fitCoeffs(4)*fitCoeffs(3), sqrt((fitCoeffs(1)*d(3))^2 + (fitCoeffs(3)*d(1))^2 + d(2)^2 + (fitCoeffs(4)*d(3))^2 + (fitCoeffs(3)*d(4))^2));
end

if (OutputFile == 1)
    for k = 1:(length(fitCoeffs)-1)
        fprintf(fout, '%f, ', fitCoeffs(k));
        fprintf(fout, '%f, ', d(k));
    end
    fprintf(fout, '%f, ', fitCoeffs(length(fitCoeffs)));
    fprintf(fout, '%f\n', d(length(fitCoeffs)));
end

%% Calculate and Plot Constraints

% Calculate using functions defined at end of script
phaseConstraints = calc_dfdot_phase_constraint(10.^imat, (10.^Tobsmat).*3600, fdot_sig);
freqConstraints = calc_dfdot_freq_constraint(10.^imat, (10.^Tobsmat).*3600, fdot_sig);

% Add to figure 1000
figure(1000);
hold on;
sf10001 = surf(Tobsmat, imat, log10(phaseConstraints));
sf10002 = surf(Tobsmat, imat, log10(freqConstraints));
sf10001.EdgeAlpha = 0.25;
sf10002.EdgeAlpha = 0.25;
sf10001.FaceAlpha = 0.25;
sf10002.FaceAlpha = 0.25;
sf10001.FaceColor = 'blue';
sf10002.FaceColor = 'red';
if (fitDegree == pi || fitDegree == 2*pi)
    legend('Calculated Data', 'Fit', 'Intersection', 'Phase Constraint', 'Frequency Constraint');
else
    legend('Calculated Data', 'Fit', 'Phase Constraint', 'Frequency Constraint');
end
hold off;

for k = 1:length(Tobsvec)
    fprintf('Range of template spacings for i = %d: %e\n', Nsegvec(k), range(data(k:end, Nsegvec(k))))
end

%% Functions

function phaseConstraint = calc_dfdot_phase_constraint(i, Tobs, fdot_sig, desiredPhaseMismatch)
    if (nargin == 3)
        desiredPhaseMismatch = 0.20;
    end
    Tcoh = sqrt(0.5./abs(fdot_sig));
    Tcoh_eff = i.*floor(Tcoh);
    Nseg_eff = Tobs./Tcoh_eff;
    phaseConstraint = 3.*sqrt(5.*desiredPhaseMismatch)./(pi.*Tcoh_eff.^2)./sqrt(Nseg_eff);
end

function freqConstraint = calc_dfdot_freq_constraint(i, Tobs, fdot_sig, desiredFreqMismatch)
    if (nargin == 3)
        desiredFreqMismatch = 1;
    end
    Tcoh = sqrt(0.5./abs(fdot_sig));
    binWidth_eff = (i.*Tcoh).^(-1);
    freqConstraint = desiredFreqMismatch*binWidth_eff./Tobs;
end

