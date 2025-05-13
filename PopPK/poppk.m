clear; clc;

% --- Load population data ---
load('valsartan_random_pop_case1_case2.mat');  % weights, receptors_case1, receptors_case2

%% --- Case 1: healthy baseline receptor concentration ---
n = length(weights);
AUC_angII_AT1R_case1 = zeros(n,1);  % output array

for i = 1:n
    % Set individual parameters
    p.Vd = weights(i) * 0.6;                 % L
    p.k_CL_1 = 2 / p.Vd / 0.05;              % 1/h
    p.k_a = 1.409;                           % 1/h
    p.C0_2 = receptors_case1(i);             % receptor concentration (µM)
    p.C0_3 = 0.0239;                         % baseline AngII concentration (µM)
    p.doses = [0, 80];                       % single 80mg dose
    p.a = 1;                                 % feedback on

    % Run simulation
    try
        [t, y] = sim0_v2(p);
        AUC_angII_AT1R_case1(i) = trapz(t, y(:,5));  % AUC of AngII–R complex
    catch
        AUC_angII_AT1R_case1(i) = NaN;
        warning('Simulation failed for subject %d', i);
    end
end

%% --- Case 2: Disease baseline receptor concentration ---
n = length(weights);
AUC_angII_AT1R_case2 = zeros(n,1);  % output array

for i = 1:n
    % Set individual parameters
    p.Vd = weights(i) * 0.6;                 % L
    p.k_CL_1 = 2 / p.Vd / 0.05;              % 1/h
    p.k_a = 1.409;                           % 1/h
    p.C0_2 = receptors_case2(i);             % receptor concentration (µM)
    p.C0_3 = 0.0239;                         % baseline AngII concentration (µM)
    p.doses = [0, 80];                       % single 80mg dose
    p.a = 1;                                 % feedback on

    % Run simulation
    try
        [t, y] = sim0_v2(p);
        AUC_angII_AT1R_case2(i) = trapz(t, y(:,5));  % AUC of AngII–R complex
    catch
        AUC_angII_AT1R_case2(i) = NaN;
        warning('Simulation failed for subject %d', i);
    end
end

% --- Determine uniform color and y-axis limits ---
min_receptor = min([receptors_case1; receptors_case2]);
max_receptor = max([receptors_case1; receptors_case2]);

min_auc = min([AUC_angII_AT1R_case1; AUC_angII_AT1R_case2]);
max_auc = max([AUC_angII_AT1R_case1; AUC_angII_AT1R_case2]);

% % --- Scatter plot with receptor concentration as color ---
% figure;
% scatter(weights, AUC_angII_AT1R_case1, 60, receptors_case1, 'filled');
% xlabel('Weight (kg)');
% ylabel('AUC of AngII–AT1R Complex (\muM·h)');
% title('Weight & Baseline [AT1R] vs AUC (Healthy)');
% colormap('parula');
% clim([min_receptor, max_receptor]);  % uniform color scale
% ylim([min_auc-1e-4, max_auc+1e-4]);             % uniform y-axis
% c = colorbar;
% c.Label.String = 'Receptor Concentration (\muM)';
% grid on;
% saveas(gcf, 'Case1_Weight_Receptor_Vary.png');
% 
% % --- Scatter plot with receptor concentration as color ---
% figure;
% scatter(weights, AUC_angII_AT1R_case2, 60, receptors_case2, 'filled');
% xlabel('Weight (kg)');
% ylabel('AUC of AngII–AT1R Complex (\muM·h)');
% title('Weight & Baseline [AT1R] vs AUC (Diseased)');
% colormap('parula');
% clim([min_receptor, max_receptor]);  % same as Case 1
% ylim([min_auc-1e-4, max_auc+1e-4]);             % same as Case 1
% c = colorbar;
% c.Label.String = 'Receptor Concentration (\muM)';
% grid on;
% saveas(gcf, 'Case2_Diseased_Receptor_High.png');
% 
% % Save results
% save('valsartan_auc_angII_case2.mat', 'weights', 'receptors_case2', 'AUC_angII_AT1R_case2');
% disp('AUC results for Case 2 saved.');

%% --- Case 3: Vary weight only, fixed receptor ---
AUC_case3 = zeros(n,1);
fixed_receptor = 0.08585;

for i = 1:n
    p.Vd = weights(i) * 0.6;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.k_a = 1.409;
    p.C0_2 = fixed_receptor;
    p.C0_3 = 0.0239;
    p.doses = [0, 80];
    p.a = 1;

    try
        [t, y] = sim0_v2(p);
        AUC_case3(i) = trapz(t, y(:,5));
    catch
        AUC_case3(i) = NaN;
    end
end

%% --- Case 4: Fixed weight, vary receptor only ---
AUC_case4 = zeros(n,1);
fixed_weight = 73.9;
fixed_Vd = fixed_weight * 0.6;

for i = 1:n
    p.Vd = fixed_Vd;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.k_a = 1.409;
    p.C0_2 = receptors_case1(i);  % Varying only receptor
    p.C0_3 = 0.0239;
    p.doses = [0, 80];
    p.a = 1;

    try
        [t, y] = sim0_v2(p);
        AUC_case4(i) = trapz(t, y(:,5));
    catch
        AUC_case4(i) = NaN;
    end
end

% %% --- Plotting: Case 3 (Weight only) ---
% figure;
% scatter(weights, AUC_case3, 60, 'filled');
% xlabel('Weight (kg)');
% ylabel('AUC of AngII–AT1R Complex (\muM·h)');
% title('Case 3: Weight Variation Only (Fixed Receptor)');
% grid on;
% saveas(gcf, 'Case3_Weight_Only.png');
% 
% %% --- Plotting: Case 4 (Receptor only) ---
% figure;
% scatter(receptors_case1, AUC_case4, 60, 'filled');
% xlabel('Receptor Concentration (\muM)');
% ylabel('AUC of AngII–AT1R Complex (\muM·h)');
% title('Case 4: Receptor Variation Only (Fixed Weight)');
% grid on;
% saveas(gcf, 'Case4_Receptor_Only.png');

%% Save results
save('valsartan_auc_angII_case3_case4.mat', 'weights', 'receptors_case1', 'AUC_case3', 'AUC_case4');
disp('AUC results for Case 3 and 4 saved.');

%% --- Case 5: Vary weight only, fixed receptor based on disease population mean ---
AUC_case5 = zeros(n,1);
fixed_receptor = 0.08585 * 1.9;

for i = 1:n
    p.Vd = weights(i) * 0.6;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.k_a = 1.409;
    p.C0_2 = fixed_receptor;
    p.C0_3 = 0.0239;
    p.doses = [0, 80];
    p.a = 1;

    try
        [t, y] = sim0_v2(p);
        AUC_case5(i) = trapz(t, y(:,5));
    catch
        AUC_case5(i) = NaN;
    end
end

%% --- Case 6: Fixed weight, vary receptor only, disease population ---
AUC_case6 = zeros(n,1);
fixed_weight = 73.9;
fixed_Vd = fixed_weight * 0.6;

for i = 1:n
    p.Vd = fixed_Vd;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.k_a = 1.409;
    p.C0_2 = receptors_case2(i);  % Varying only receptor
    p.C0_3 = 0.0239;
    p.doses = [0, 80];
    p.a = 1;

    try
        [t, y] = sim0_v2(p);
        AUC_case6(i) = trapz(t, y(:,5));
    catch
        AUC_case6(i) = NaN;
    end
end

%% --- Plotting: Case 3 and 5 (Weight only, healthy and diseased population) ---
figure;
hold on;

% Plot Case 3: Healthy population receptor (fixed)
scatter(weights, AUC_case3, 60, 'o', 'MarkerEdgeColor',[0.4, 0.7, 1], 'MarkerFaceColor', 'none');

% Plot Case 5: Diseased population receptor (fixed)
scatter(weights, AUC_case5, 60, 'o', 'MarkerEdgeColor',[1, 0.5, 0.5], 'MarkerFaceColor', 'none');

xlabel('Weight (kg)', 'FontSize', 15);
ylabel('AUC of AngII–AT1R Complex (\muM·h)', 'FontSize', 15);
%title('Weight vs AUC (Fixed baseline [AT1R] across the same group)', 'FontSize', 15);
legend('Healthy', 'Diseased', 'Location', 'southeast', 'FontSize', 15);
ylim([1e-5, 12e-4]);
set(gca, 'FontSize', 15); % Set axes tick labels to 15pt
grid on;

saveas(gcf, 'Varied_Weight_AUC.png');

%% --- Plotting: Case 4 and 6 (Receptor only, healthy and diseased population) ---
figure;
hold on;

% Plot Case 4: Healthy receptor distribution
scatter(receptors_case1, AUC_case4, 60, 'o', 'MarkerEdgeColor',[0.4, 0.7, 1], 'MarkerFaceColor', 'none');

% Plot Case 6: Diseased receptor distribution
scatter(receptors_case2, AUC_case6, 60, 'o', 'MarkerEdgeColor',[1, 0.5, 0.5], 'MarkerFaceColor', 'none');

xlabel('Receptor Concentration (\muM)', 'FontSize', 15);
ylabel('AUC of AngII–AT1R Complex (\muM·h)', 'FontSize', 15);
%title('Receptor vs AUC (Fixed weight)', 'FontSize', 15);
legend('Healthy', 'Diseased', 'Location', 'southeast', 'FontSize', 15);
ylim([1e-5, 12e-4]);
set(gca, 'FontSize', 15); % Set axes tick labels to 15pt
grid on;

saveas(gcf, 'Varied_Receptor_AUC.png');


%% Save results
save('valsartan_auc_angII_case5_case6.mat', 'weights', 'receptors_case2', 'AUC_case5', 'AUC_case6');
disp('AUC results for Case 5 and 6 saved.');
