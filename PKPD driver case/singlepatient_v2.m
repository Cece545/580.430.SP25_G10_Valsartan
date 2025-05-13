close all;
clear;

%% Step 1: Set parameters
p.C0_2 = 0.08585;     % Initial Receptor concentration
p.C0_3 = 0.0239;      % Initial Angiotensin II concentration
% Dosing every 24 hours for 8 doses of 80 mg each
num_doses = 1;
dose_interval = 24; % hours
p.dose_amount = 80;   % mg

p.doses = zeros(num_doses, 2); % Preallocate [time, dose] matrix

for i = 1:num_doses
    p.doses(i, 1) = (i-1) * dose_interval; % Time of dose (in hours)
    p.doses(i, 2) = dose_amount;           % Dose amount (in mg)
end

%% Step 2: Run simulation
[time, y] = sim0_v2(p);  % y: matrix of size [time x 9]

%% Step 3: Variable labels (for clarity in plots)
labels = {
    '1. Free Valsartan', ...
    '2. Free Receptors', ...
    '3. Angiotensin II', ...
    '4. Valsartan-Receptor Complex', ...
    '5. Ang II-Receptor Complex', ...
    '6. Valsartan-Albumin Complex', ...
    '7. Cleared Valsartan', ...
    '8. Cleared Angiotensin II', ...
    '9. Ang II Production (Virtual)', ...
    '10. Valsartan in Gut'};

%% Step 4: Plot all 9 variables
figure('Name','All State Variables Over Time','NumberTitle','off');
for i = 1:10
    subplot(3,4,i);  % 3x4 grid of subplots
    plot(time, y(:,i), 'LineWidth', 2);
    title(labels{i});
    xlabel('Time (hours)');
    ylabel('Concentration or Amount');
    grid on;
end

sgtitle('Simulation Output: All State Variables');
