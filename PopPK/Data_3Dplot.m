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

% Save Case 1 results to CSV
case1_data = table(weights, receptors_case1, AUC_angII_AT1R_case1, 'VariableNames', {'Weight_kg', 'Receptor_Healthy_uM', 'AUC_AngII_AT1R_Healthy'});
writetable(case1_data, 'valsartan_case1_healthy.csv');
disp('Case 1 results saved to valsartan_case1_healthy.csv');

%% --- Case 2: Disease baseline receptor concentration ---
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

% Save Case 2 results to CSV
case2_data = table(weights, receptors_case2, AUC_angII_AT1R_case2, 'VariableNames', {'Weight_kg', 'Receptor_Diseased_uM', 'AUC_AngII_AT1R_Diseased'});
writetable(case2_data, 'valsartan_case2_diseased.csv');
disp('Case 2 results saved to valsartan_case2_diseased.csv');
