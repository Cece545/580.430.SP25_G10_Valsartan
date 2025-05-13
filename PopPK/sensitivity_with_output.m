clear; close all; clc;

%% Load Simulated Population
load('valsartan_random_pop_case1_case2.mat');  % Contains weights and receptors_case1
n_patients = length(weights);

%% Initialize relative sensitivity arrays
S_Cmax_weight = zeros(n_patients,1);
S_Cmin_weight = zeros(n_patients,1);
S_AUC_weight  = zeros(n_patients,1);

S_Cmax_receptor = zeros(n_patients,1);
S_Cmin_receptor = zeros(n_patients,1);
S_AUC_receptor  = zeros(n_patients,1);

%% Constants
base_dose = 80;
bioavailability = 0.25;
Vd_per_kg = 0.6;
k_a = 1.409;
C0_3 = 0.0239;

for i = 1:n_patients
    %% Baseline parameters
    p = struct();
    p.C0_2 = receptors_case1(i);
    p.C0_3 = C0_3;
    p.weight = weights(i);
    p.num_doses = 1;
    p.dose_interval = 24;
    p.dose_amount = base_dose;
    p.doses = [0, base_dose];
    p.k_a = k_a;
    p.Vd = p.weight * Vd_per_kg;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.production = 0;
    p.a = 1;

    %% Baseline simulation
    [time_base, y_base] = sim0_v2(p);
    valsartan_base = y_base(:,1);
    Cmax_base = max(valsartan_base);
    mask_last24 = time_base >= max(time_base) - 24;
    Cmin_base = mean(valsartan_base(mask_last24));
    AUC_base = trapz(time_base, valsartan_base);

    %% Perturb weight
    delta_w = 0.05 * p.weight;
    p_w = p;
    p_w.weight = p.weight + delta_w;
    p_w.Vd = p_w.weight * Vd_per_kg;
    p_w.k_CL_1 = 2 / p_w.Vd / 0.05;
    [time_w, y_w] = sim0_v2(p_w);
    valsartan_w = y_w(:,1);
    Cmax_w = max(valsartan_w);
    Cmin_w = mean(valsartan_w(time_w >= max(time_w) - 48));
    AUC_w = trapz(time_w, valsartan_w);

    % Relative sensitivity (dimensionless)
    S_Cmax_weight(i) = ((Cmax_w - Cmax_base) / Cmax_base) / (delta_w / p.weight);
    S_Cmin_weight(i) = ((Cmin_w - Cmin_base) / Cmin_base) / (delta_w / p.weight);
    S_AUC_weight(i)  = ((AUC_w  - AUC_base)  / AUC_base)  / (delta_w / p.weight);

    %% Perturb receptor concentration
    delta_r = 0.05 * p.C0_2;
    p_r = p;
    p_r.C0_2 = p.C0_2 + delta_r;
    [time_r, y_r] = sim0_v2(p_r);
    valsartan_r = y_r(:,1);
    Cmax_r = max(valsartan_r);
    Cmin_r = mean(valsartan_r(time_r >= max(time_r) - 24));
    AUC_r = trapz(time_r, valsartan_r);

    % Relative sensitivity (dimensionless)
    S_Cmax_receptor(i) = ((Cmax_r - Cmax_base) / Cmax_base) / (delta_r / p.C0_2);
    S_Cmin_receptor(i) = ((Cmin_r - Cmin_base) / Cmin_base) / (delta_r / p.C0_2);
    S_AUC_receptor(i)  = ((AUC_r  - AUC_base)  / AUC_base)  / (delta_r / p.C0_2);
end

%% === HEATMAP for First 20 Patients ===
n_show = 100;
heatmap_data = [ ...
    S_Cmax_weight(1:n_show), ...
    S_Cmin_weight(1:n_show), ...
    S_AUC_weight(1:n_show), ...
    S_Cmax_receptor(1:n_show), ...
    S_Cmin_receptor(1:n_show), ...
    S_AUC_receptor(1:n_show)];

metric_labels = {'Rel dC_{max}/dWeight', 'Rel dC_{min}/dWeight', 'Rel dAUC/dWeight', ...
                 'Rel dC_{max}/dReceptor', 'Rel dC_{min}/dReceptor', 'Rel dAUC/dReceptor'};
patient_labels = arrayfun(@(x) sprintf('Patient %d', x), 1:n_show, 'UniformOutput', false);

figure('Name','Relative Sensitivity Heatmap','NumberTitle','off');
h = heatmap(metric_labels, patient_labels, heatmap_data);
h.Title = 'Relative Sensitivity of Valsartan PK Metrics';
h.XLabel = 'Sensitivity Metric';
h.YLabel = 'Simulated Patients';

% Red–white–blue diverging colormap
n = 256;
half = round(n/2);
cmap = zeros(n,3);
cmap(1:half,1) = linspace(0,1,half);   % Red increases
cmap(1:half,2) = linspace(0,1,half);   % Green increases
cmap(1:half,3) = 1;                    % Blue max
cmap(half+1:end,1) = 1;                % Red max
cmap(half+1:end,2) = linspace(1,0,n-half);
cmap(half+1:end,3) = linspace(1,0,n-half);
colormap(cmap);

% Symmetric color limits centered at zero
clim = max(abs(heatmap_data(:)));
h.ColorLimits = [-clim, clim];

%% === BAR PLOT of Mean Relative Sensitivity ===
mean_sensitivities = [ ...
    mean(S_Cmax_weight), ...
    mean(S_Cmin_weight), ...
    mean(S_AUC_weight), ...
    mean(S_Cmax_receptor), ...
    mean(S_Cmin_receptor), ...
    mean(S_AUC_receptor)];

figure('Name','Mean Relative Sensitivity','NumberTitle','off');
bar(mean_sensitivities);
set(gca, 'XTickLabel', metric_labels, 'XTickLabelRotation', 45);
ylabel('Mean Relative Sensitivity (dimensionless)');
title('Average Relative Sensitivity of PK Metrics');
grid on;
save('valsartan_relative_sensitivity.mat', ...
     'S_Cmax_weight', 'S_Cmin_weight', 'S_AUC_weight', ...
     'S_Cmax_receptor', 'S_Cmin_receptor', 'S_AUC_receptor');
