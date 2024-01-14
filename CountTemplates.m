% Calculates template count of semi-coherent search over f and fdot for various i and Tobs

clear;
clc;
close all;

% Choose Settings
PlotValidiOnly = 1;
LogPlot = 1;

%% Select rectangular region of parameter space to search

% Choose min and max frequencies (in Hz)
f1 = 10;
f2 = 1000;

% Choose min and max frequency derivatives (Hz/s)
% NOTE: fdot1 < fdot2, but |fdot1| > |fdot2|, as fdot < 0
fdot1 = -5.e-7;
fdot2 = -5.e-9;

%% Calculate template count for one chosen i and Tobs

% Choose grouping number i and observation period (hr) for calculation of
% single template count
i = 100;
Tobs_hr = 32;

% Load parameters and calculate template count from analytically integated formula
params = readmatrix('SNRgModelParams.csv');

N_T = count_templates_max(i, Tobs_hr, f1, f2, fdot1, fdot2, params);

% Print template count
fprintf('For searching the space %1.f Hz <= f<= %1.f Hz \nand %1.e Hz/s <= fdot <= %1.e Hz/s\n', f1, f2, fdot1, fdot2);
fprintf('with Tobs = %1.f hr and i = %d,\n', Tobs_hr, i);
fprintf('%e templates are required.\n', N_T);
fprintf('Template Count per unit frequency: %e\n', N_T/(f2-f1));

%% Calculate Template Count for range of i and Tobs for valid constant i based on fdot2

% Create grids of Tobs and i
Tobsvec_hr = 4:0.4:40;
ivec = 1:round(max(Tobsvec_hr)*3600*sqrt(2*abs(-5.e-5)));

[Tobsmat, imat] = meshgrid(Tobsvec_hr, ivec);

% Count templates for each combination
N_Tmat = count_templates(imat, Tobsmat, f1, f2, fdot1, fdot2, params);
if (LogPlot == 1)
    N_Tmat = log10(N_Tmat);
end

% Create logical arrays to index where template counts are valid
% Based on if i exists for given Tobs and fdot
valid = is_valid_i(imat, Tobsmat, fdot2);

% Separate valid and invalid i points
N_Tmat_valid = N_Tmat;
N_Tmat_valid(~valid) = NaN;

N_Tmat_invalid = N_Tmat;
N_Tmat_invalid(valid) = NaN;


% Plot on Tobs-i plane
red = zeros(height(N_Tmat), width(N_Tmat), 3);
red(:, :, 1) = 0.75;

figure;
sf1 = surf(Tobsmat, imat, N_Tmat_valid, 'EdgeAlpha', 0);
hold on;
if (PlotValidiOnly == 0)
sf2 = surf(Tobsmat, imat, N_Tmat_invalid, red, 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
legend('Valid i', 'Invalid i');
end
xlabel('Tobs (hr)');
ylabel('i');
if (LogPlot == 1)
    zlabel('Base 10 log of Template Count');
else
    zlabel('Template Count');
end
title(sprintf('Template Counts for f on [%1.f, %1.f] Hz and fdot on [%1.e, %1.e] Hz/s', f1, f2, fdot1, fdot2));
ax = gca;
ax.FontSize = 15;
ax.LineWidth = 3;
grid on;
hold off;

%% Calculate Template Counts for range of i and Tobs using greatest possible i for each fdot range up to imax
N_Tmat_max = count_templates_max(imat, Tobsmat, f1, f2, fdot1, fdot2, params);

%% Plot 
figure;
sf3 = surf(Tobsmat, imat, N_Tmat_max, 'EdgeAlpha', 0);
zlabel('Template Count $\tilde{N}_T$', 'Interpreter', 'LaTex');
xlabel('$T_{obs}$ (hr)', 'Interpreter', 'LaTex');
ylabel('$g_{max}$', 'Interpreter', 'LaTex');
title(sprintf('Template Counts for %s on [%1.f, %1.f] Hz and %s on [%1.e, %1.e] Hz/s \nusing maximum %s up to %s', '$f_0$', f1, f2, '$\dot{f}$', fdot1, fdot2, '$g$', '$g_{max}$'), 'Interpreter', 'LaTex');
%colorbar;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
grid on;
hold off;

%% Functions (Tobs in hr for all)

% Determines if i is defined for given Tobs and fdot2
% NOT VECTORIZED
function validity = is_valid_i_scalar(i, Tobs, fdot2)
    epsilon = 1.e-12;
    i_raw = Tobs*3600*sqrt(2*abs(fdot2));
    if ((ceil(i_raw) - i_raw) <= epsilon)
        validity = i <= ceil(i_raw);
    else
        validity = i <= floor(i_raw);
    end
end

% Determines if i is defined for given Tobs and fdot2
% i and Tobs must be same size
% VECTORIZED
function validity = is_valid_i(i, Tobs, fdot2)
    validity = false(size(i));
    for k = 1:numel(validity)
        validity(k) = is_valid_i_scalar(i(k), Tobs(k), fdot2);
    end
end

% Counts templates for search over f-fdot range with given i and Tobs_in
% i_in and Tobs must be same size
function N_T = count_templates(i_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)       
        N_T = 10.^(-params(1, 4)-params(2, 4)) .* ...
              Tobs_in.^(-params(1, 1)-params(2, 1)) .* ...
              i_in.^(-params(1, 2)-params(2, 2)) .* ...
              (f2_in - f1_in) .* ...
              (abs(fdot1_in).^(1-params(1, 3)-params(2, 3)) - abs(fdot2_in).^(1-params(1, 3)-params(2, 3))) ./ ...
              (1-params(1, 3)-params(2, 3));
end

% Counts templates using max possible i for each subrange of fdot
% NOT VECTORIZED
function N_T = count_templates_max_scalar(imax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    fdot_max = -imax_in^2/2/(3600*Tobs_in)^2;
    if (is_valid_i_scalar(imax_in, Tobs_in, fdot2_in))
        N_T = count_templates(imax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    elseif (fdot_max < fdot1_in)
        N_T = count_templates_max_scalar(imax_in-1, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    else
        N_T = count_templates_max_scalar(imax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot_max, params) + ... 
              count_templates_max_scalar(imax_in-1, Tobs_in, f1_in, f2_in, fdot_max, fdot2_in, params);
    end
end

% Counts templates using max possible i for each subrange of fdot
% i_in and Tobs must be same size
% VECTORIZED
function N_T = count_templates_max(imax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    N_T = zeros(size(imax_in));
    for k = 1:numel(imax_in)
        N_T(k) = count_templates_max_scalar(imax_in(k), Tobs_in(k), f1_in, f2_in, fdot1_in, fdot2_in, params);
    end
end
