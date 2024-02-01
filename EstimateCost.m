% Calculates template count of semi-coherent search over f and fdot for various i and Tobs

clear;
clc;
close all;

% Choose Settings
PlotValidiOnly = 1;
LogPlot = 1;
ColorMapView = 1;

%% Select rectangular region of parameter space to search

% Choose min and max frequencies (in Hz)
f1 = 10;
f2 = 10.005;

% Choose min and max frequency derivatives (Hz/s)
% NOTE: fdot1 < fdot2, but |fdot1| > |fdot2|, as fdot < 0
fdot1 = -1.e-8;
fdot2 = -5.e-9;

%% Calculate template count for one chosen i and Tobs

% Choose grouping number i and observation period (hr) for calculation of
% single template count
g = 5;
Tobs_hr = 10;

% Load parameters and calculate template count from analytically integated formula
params = readmatrix('SNRgModelParams.csv');

N_T = count_templates_max(g, Tobs_hr, f1, f2, fdot1, fdot2, params);
T_comp = time_search_max(g, Tobs_hr, f1, f2, fdot1, fdot2, params);

% Print template count
fprintf('For searching the space %1.f Hz <= f<= %1.f Hz \nand %1.e Hz/s <= fdot <= %1.e Hz/s\n', f1, f2, fdot1, fdot2);
fprintf('with Tobs = %1.f hr and g = %d,\n', Tobs_hr, g);
fprintf('%e templates are required.\n', N_T);
fprintf('Template Count per unit frequency: %e (1/Hz)\n', N_T/(f2-f1));
fprintf('The computation time is %e s\n', T_comp);
fprintf('Time per unit frequency: %e (s/Hz)\n', T_comp/(f2-f1));

%% Calculate Template Count for range of i and Tobs for valid constant i based on fdot2

% Create grids of Tobs and g
Tobsvec_hr = 1:0.4:40;
Nsegline = floor(3600.*Tobsvec_hr.*sqrt(2.*abs(fdot1)));
gvec = 1:round(max(Tobsvec_hr)*3600*sqrt(2*abs(-5.e-5)));

[Tobsmat, gmat] = meshgrid(Tobsvec_hr, gvec);

% Count templates for each combination
N_Tmat = count_templates(gmat, Tobsmat, f1, f2, fdot1, fdot2, params);
if (LogPlot == 1)
    N_Tmat = log10(N_Tmat);
end

% Create logical arrays to index where template counts are valid
% Based on if g exists for given Tobs and fdot
valid = is_valid_g(gmat, Tobsmat, fdot2);

% Separate valid and invalid i points
N_Tmat_valid = N_Tmat;
N_Tmat_valid(~valid) = NaN;

N_Tmat_invalid = N_Tmat;
N_Tmat_invalid(valid) = NaN;


% Plot on Tobs-i plane
red = zeros(height(N_Tmat), width(N_Tmat), 3);
red(:, :, 1) = 0.75;

figure;
sf1 = surf(Tobsmat, gmat, N_Tmat_valid, 'EdgeAlpha', 0);
hold on;
if (PlotValidiOnly == 0)
sf2 = surf(Tobsmat, gmat, N_Tmat_invalid, red, 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
legend('Valid g', 'Invalid g');
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

model_funcng = @(a, fdot_log, Tobs_log) a(1).*fdot_log + ...
                                        a(2)*Tobs_log + a(3);

fprintf('Time for single template: %e s\n', 10.^model_funcng(params(4, :), log10(5.e-5), log10(40)));

%% Calculate Template Counts for range of g and Tobs using greatest possible g for each fdot range up to gmax
N_Tmat_max = count_templates_max(gmat, Tobsmat, f1, f2, fdot1, fdot2, params);

%% Plot 
figure;
hold on;
sf3 = surf(Tobsmat, gmat, N_Tmat_max, 'EdgeAlpha', 0);
zlabel('Template Count $\tilde{N}_T$', 'Interpreter', 'LaTex');
xlabel('$T_{obs}$ (hr)', 'Interpreter', 'LaTex');
ylabel('$g_{max}$', 'Interpreter', 'LaTex');
title(sprintf('Template Counts for %s on [%1.f, %1.f] Hz and %s on [%1.e, %1.e] Hz/s \nusing maximum %s up to %s', '$f_0$', f1, f2, '$\dot{f}$', fdot1, fdot2, '$g$', '$g_{max}$'), 'Interpreter', 'LaTex');
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
ax.ZScale = 'log';
set(gca,'ColorScale','log')
if (ColorMapView == 1)
    colorbar;
    p3 = plot3(Tobsvec_hr, Nsegline, max(N_Tmat_max(:)).*ones(size(Tobsvec_hr)), '-r', 'LineWidth', 3);
    legend('template counts', 'maximum non-redundant $g_{max}$', 'Interpreter', 'LaTex', 'Location', 'northwest');
    view([0, 90]);
end
grid on;
hold off;

%% Calculate Search Time for range of g and Tobs using greatest possible g for each fdot range up to gmax
T_compmat_max = time_search_max(gmat, Tobsmat, f1, f2, fdot1, fdot2, params);

%% Plot 
figure;
hold on;
sf4 = surf(Tobsmat, gmat, T_compmat_max, 'EdgeAlpha', 0);
zlabel('Computation Time $\tilde{T}_{comp} (s)$', 'Interpreter', 'LaTex');
xlabel('$T_{obs}$ (hr)', 'Interpreter', 'LaTex');
ylabel('$g_{max}$', 'Interpreter', 'LaTex');
title(sprintf('Computation Time of Search with %s on [%1.f, %1.f] Hz and %s on [%1.e, %1.e] Hz/s \nusing maximum %s up to %s', '$f_0$', f1, f2, '$\dot{f}$', fdot1, fdot2, '$g$', '$g_{max}$'), 'Interpreter', 'LaTex');
colorbar;
ax = gca;
ax.FontSize = 20;
ax.LineWidth = 3;
ax.ZScale = 'log';
set(gca,'ColorScale','log')
if (ColorMapView == 1)
    colorbar;
    p4 = plot3(Tobsvec_hr, Nsegline, max(T_compmat_max(:)).*ones(size(Tobsvec_hr)), '-r', 'LineWidth', 3);
    legend('time (s)', 'maximum non-redundant $g_{max}$', 'Interpreter', 'LaTex', 'Location', 'northwest');
    view([0, 90]);
end
grid on;
hold off;

%% Calculate Search Time for range of g and Tobs using every valid g which is a power of 2

model_func = @(a, fdot_log, Tobs_log, g_log) a(1).*Tobs_log + ...
                                              a(2).*g_log + ...
                                              a(3)*fdot_log + a(4);

model_funcX = @(a, fdot_log, Tobs_log)  a(1).*fdot_log.^2 + ...
                                        a(2).*fdot_log.*Tobs_log + ...
                                        a(3).*Tobs_log.^2 + ...
                                        a(4).*fdot_log + ...
                                        a(5).*Tobs_log + ...
                                        a(6);

rho_comp = @(fdot, Tobs_hr) 10.^(model_funcX(params(5, :), log10(abs(fdot)), log10(Tobs_hr)))./10.^(model_func(params(1, :), log10(abs(fdot)), log10(Tobs_hr), log10(2.^floor(log2(floor(3600*Tobs_hr*sqrt(2*abs(fdot))))))))./10.^(model_func(params(2, :), log10(abs(fdot)), log10(Tobs_hr), log10(2.^floor(log2(floor(3600*Tobs_hr*sqrt(2*abs(fdot))))))));
T_compmat2 = zeros(size(Tobsvec_hr));

for i = 1:length(Tobsvec_hr)
    T_compmat2(i) = (f2-f1)*integral(@(x) rho_comp(x, Tobsvec_hr(i)), fdot1, fdot2);
end

figure;plot(Tobsvec_hr, T_compmat2, 'LineWidth', 3); set(gca, 'YScale', 'log');
title('Computation Times using every valid g that is a power of 2'); xlabel('Tobs (hr)'); ylabel('Time (s)');
ax = gca; ax.FontSize = 15;

%% Functions (Tobs in hr for all)

% Determines if g is defined for given Tobs and fdot2
% NOT VECTORIZED
function validity = is_valid_g_scalar(g, Tobs, fdot2)
    epsilon = 1.e-12;
    g_raw = Tobs*3600*sqrt(2*abs(fdot2));
    if ((ceil(g_raw) - g_raw) <= epsilon)
        validity = g <= ceil(g_raw);
    else
        validity = g <= floor(g_raw);
    end
end

% Determines if g is defined for given Tobs and fdot2
% g and Tobs must be same size
% VECTORIZED
function validity = is_valid_g(g, Tobs, fdot2)
    validity = false(size(g));
    for k = 1:numel(validity)
        validity(k) = is_valid_g_scalar(g(k), Tobs(k), fdot2);
    end
end

% Counts templates for search over f-fdot range with given g and Tobs_in
% g_in and Tobs must be same size
function N_T = count_templates(g_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)       
        N_T = 10.^(-params(1, 4)-params(2, 4)) .* ...
              Tobs_in.^(-params(1, 1)-params(2, 1)) .* ...
              g_in.^(-params(1, 2)-params(2, 2)) .* ...
              (f2_in - f1_in) .* ...
              (abs(fdot1_in).^(1-params(1, 3)-params(2, 3)) - abs(fdot2_in).^(1-params(1, 3)-params(2, 3))) ./ ...
              (1-params(1, 3)-params(2, 3));
end

% Counts templates using max possible g for each subrange of fdot
% NOT VECTORIZED
function N_T = count_templates_max_scalar(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    fdot_max = -gmax_in^2/2/(3600*Tobs_in)^2;
    if (is_valid_g_scalar(gmax_in, Tobs_in, fdot2_in))
        N_T = count_templates(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    elseif (fdot_max < fdot1_in)
        N_T = count_templates_max_scalar(gmax_in-1, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    else
        N_T = count_templates_max_scalar(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot_max, params) + ... 
              count_templates_max_scalar(gmax_in-1, Tobs_in, f1_in, f2_in, fdot_max, fdot2_in, params);
    end
end

% Counts templates using max possible g for each subrange of fdot
% g_in and Tobs must be same size
% VECTORIZED
function N_T = count_templates_max(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    N_T = zeros(size(gmax_in));
    for k = 1:numel(gmax_in)
        N_T(k) = count_templates_max_scalar(gmax_in(k), Tobs_in(k), f1_in, f2_in, fdot1_in, fdot2_in, params);
    end
end

% Estimate cost for search over f-fdot range with given g and Tobs_in
% g_in and Tobs must be same size
function T_T = time_search(g_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)       
        T_T = 10.^(params(4, 3)-params(1, 4)-params(2, 4)) .* ...
              Tobs_in.^(params(4, 2)-params(1, 1)-params(2, 1)) .* ...
              g_in.^(-params(1, 2)-params(2, 2)) .* ...
              (f2_in - f1_in) .* ...
              (abs(fdot1_in).^(1+params(4, 1)-params(1, 3)-params(2, 3)) - abs(fdot2_in).^(1+params(4, 1)-params(1, 3)-params(2, 3))) ./ ...
              (1+params(4, 1)-params(1, 3)-params(2, 3));
end

% Time search using max possible g for each subrange of fdot
% NOT VECTORIZED
function T_T = time_search_max_scalar(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    fdot_max = -gmax_in^2/2/(3600*Tobs_in)^2;
    if (is_valid_g_scalar(gmax_in, Tobs_in, fdot2_in))
        T_T = time_search(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    elseif (fdot_max < fdot1_in)
        T_T = time_search_max_scalar(gmax_in-1, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params);
    else
        T_T = time_search_max_scalar(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot_max, params) + ... 
              time_search_max_scalar(gmax_in-1, Tobs_in, f1_in, f2_in, fdot_max, fdot2_in, params);
    end
end

% Time search using max possible g for each subrange of fdot
% g_in and Tobs must be same size
% VECTORIZED
function T_comp = time_search_max(gmax_in, Tobs_in, f1_in, f2_in, fdot1_in, fdot2_in, params)
    T_comp = zeros(size(gmax_in));
    for k = 1:numel(gmax_in)
        T_comp(k) = time_search_max_scalar(gmax_in(k), Tobs_in(k), f1_in, f2_in, fdot1_in, fdot2_in, params);
    end
end
