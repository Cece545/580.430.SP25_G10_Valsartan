clear; clc;

% === Load population data ===
load('valsartan_random_pop_case1_case2.mat');  % weights, receptors_case1, receptors_case2

%% --- Pre-allocate arrays ---
n = length(weights);
AUC_angII_AT1R_case1 = NaN(n,1);
AUC_angII_AT1R_case2 = NaN(n,1);

%% === Run simulations for both cases ===
for i = 1:n
    % Common parameters
    p.Vd = weights(i) * 0.6;
    p.k_CL_1 = 2 / p.Vd / 0.05;
    p.k_a = 1.409;
    p.C0_3 = 0.0239;
    p.doses = [0, 80];
    p.a = 1;

    % --- Case 1: healthy receptor concentration ---
    p.C0_2 = receptors_case1(i);
    try
        [t, y] = sim0_v2(p);
        AUC_angII_AT1R_case1(i) = trapz(t, y(:,5));
    catch
        warning('Case 1 failed for subject %d', i);
    end

    % --- Case 2: disease receptor concentration ---
    p.C0_2 = receptors_case2(i);
    try
        [t, y] = sim0_v2(p);
        AUC_angII_AT1R_case2(i) = trapz(t, y(:,5));
    catch
        warning('Case 2 failed for subject %d', i);
    end
end

%% === Combine all valid data to unify axis and colorbar ===
valid1 = ~isnan(AUC_angII_AT1R_case1);
valid2 = ~isnan(AUC_angII_AT1R_case2);
all_weights = [weights(valid1); weights(valid2)];
all_receptors = [receptors_case1(valid1); receptors_case2(valid2)];
all_AUCs = [AUC_angII_AT1R_case1(valid1); AUC_angII_AT1R_case2(valid2)];

x_range = linspace(min(all_weights), max(all_weights), 100);
y_range = linspace(min(all_receptors), max(all_receptors), 100);
c_min = min(all_AUCs);
c_max = max(all_AUCs);

%% === Case 1 Plot ===
[wq1, rq1] = meshgrid(x_range, y_range);
AUC_interp1 = griddata(weights(valid1), receptors_case1(valid1), AUC_angII_AT1R_case1(valid1), wq1, rq1, 'natural');

figure;
contourf(wq1, rq1, AUC_interp1, 20, 'LineColor', 'none');
xlabel('Weight (kg)');
ylabel('Receptor Concentration (\muM)');
title('Case 1: AUC of AngII–AT1R Complex');
colormap('parula');
c = colorbar;
c.Label.String = 'AUC (\muM·h)';
caxis([c_min c_max]);
xlim([min(x_range) max(x_range)]);
ylim([min(y_range) max(y_range)]);
grid on;

%% === Case 2 Plot ===
[wq2, rq2] = meshgrid(x_range, y_range);
AUC_interp2 = griddata(weights(valid2), receptors_case2(valid2), AUC_angII_AT1R_case2(valid2), wq2, rq2, 'natural');

figure;
contourf(wq2, rq2, AUC_interp2, 20, 'LineColor', 'none');
xlabel('Weight (kg)');
ylabel('Receptor Concentration (\muM)');
title('Case 2: AUC of AngII–AT1R Complex');
colormap('parula');
c = colorbar;
c.Label.String = 'AUC (\muM·h)';
caxis([c_min c_max]);
xlim([min(x_range) max(x_range)]);
ylim([min(y_range) max(y_range)]);
grid on;
saveas(gcf, 'Case2_Contour_Weight_Receptor_AUC.png');
%% === Scatter Plot: Case 1 ===
figure;
scatter(weights(valid1), AUC_angII_AT1R_case1(valid1), 40, receptors_case1(valid1), 'filled');
xlabel('Weight (kg)');
ylabel('AUC of AngII–AT1R Complex (\muM·h)');
title('a.');
colormap('parula');
c = colorbar;
c.Label.String = 'Receptor Concentration (\muM)';
caxis([min(all_receptors) max(all_receptors)]);
xlim([min(x_range) max(x_range)]);
ylim([min(all_AUCs) max(all_AUCs)]);
grid on;

%% === Scatter Plot: Case 2 ===
figure;
scatter(weights(valid2), AUC_angII_AT1R_case2(valid2), 40, receptors_case2(valid2), 'filled');
xlabel('Weight (kg)');
ylabel('AUC of AngII–AT1R Complex (\muM·h)');
title('b.');
colormap('parula');
c = colorbar;
c.Label.String = 'Receptor Concentration (\muM)';
caxis([min(all_receptors) max(all_receptors)]);
xlim([min(x_range) max(x_range)]);
ylim([min(all_AUCs) max(all_AUCs)]);
grid on;
