% valsartan_generate_random_population.m
% Generate 100 simulated individuals with random weight and receptor expression

rng(42);  % For reproducibility

% === Generate weights (kg) ===
n = 100;
weight_mean = 73.9;
weight_sd = 7;
weights = normrnd(weight_mean, weight_sd, [n, 1]);

% === Case 1: Healthy population receptor concentration ===
mean_r1 = 0.08585;
rel_sd_r1 = 0.24 / 0.80;
sd_r1 = mean_r1 * rel_sd_r1;a
receptors_case1 = normrnd(mean_r1, sd_r1, [n, 1]);

% === Case 2: Diseased population receptor concentration ===
mean_r2 = mean_r1 * 1.9;
rel_sd_r2 = 0.27 / 1.49;
sd_r2 = mean_r2 * rel_sd_r2;
receptors_case2 = normrnd(mean_r2, sd_r2, [n, 1]);

% === Randomly shuffle receptor data ===
receptors_case1 = receptors_case1(randperm(n));
receptors_case2 = receptors_case2(randperm(n));

% === Save to .mat ===
save('valsartan_random_pop_case1_case2.mat', 'weights', 'receptors_case1', 'receptors_case2');

disp('Dataset saved as valsartan_random_pop_case1_case2.mat');
