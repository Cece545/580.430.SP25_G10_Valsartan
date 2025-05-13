close all;
clear;

%% Step 1: Set parameters
p.C0_2 = 0.08585;     % Initial Receptor concentration
p.C0_3 = 0.0239;      % Initial Angiotensin II concentration
% Dosing every 24 hours for 8 doses of 80 mg each
p.num_doses = 1;
p.dose_interval = 24; % hours
p.dose_amount = 80;   % mg
p.doses = zeros(p.num_doses, 2); % Preallocate [time, dose] matrix


for i = 1:p.num_doses
    p.doses(i, 1) = (i-1) * p.dose_interval; % Time of dose (in hours)
    p.doses(i, 2) = p.dose_amount;           % Dose amount (in mg)
end

%% Step 2: Run simulation
p.production = 0;
[time, y] = sim0_v3(p);

%% Step 3: Variable labels (for clarity in plots)
% labels = {
%     '1. Free Valsartan', ...
%     '2. Free Receptors', ...
%     '3. Free Angiotensin II', ...
%     '4. Valsartan-Receptor Complex', ...
%     '5. Ang II-Receptor Complex', ...
%     '6. Valsartan-Albumin Complex', ...
%     '7. Cleared Valsartan', ...
%     '8. Cleared Angiotensin II', ...
%     '9. Ang II Production (Virtual)', ...
%     '10. Valsartan in Gut'};

%% Step 4: Plot variables

figure('Name','Valsartan Molecular Balance','NumberTitle','off');

hold on;

diff = (1000000/435.519 * 0.25*80-y(:,1)*17-y(:,4)*17-y(:,6)*17-y(:,7)*17-y(:,10)*17) / (1000000/435.519 * 0.25*80)

plot(time, diff,'LineWidth', 2);
hold off;

xlabel('Time (hours)');
ylabel('Mole (nmol)');
title('Valsartan Balance Over Time');
legend('w/ feedback','resting value');
grid on;

% 
% figure('Name','Ang II for Different Doses','NumberTitle','off');
% 
% hold on;
% for i = 1:length(dose_levels)
%     plot(time_storage{i}, y_storage{i}(:,3), 'LineWidth', 2);  % Column 5: Ang II-Receptor Complex
% end
% 
% hold off;
% 
% xlabel('Time (hours)');
% ylabel('Concentration');
% title('Ang II Over Time for Different Valsartan Doses');
% legend('0 mg','40 mg','80 mg','160 mg','Threshold');
% grid on;
% 
% figure('Name','Valsartan-Receptor Complex for Different Doses','NumberTitle','off');
% 
% hold on;
% for i = 1:length(dose_levels)
%     plot(time_storage{i}, y_storage{i}(:,4), 'LineWidth', 2);  % Column 5: Ang II-Receptor Complex
% end
% 
% hold off;
% 
% xlabel('Time (hours)');
% ylabel('Concentration');
% title('Valsartan-Receptor Complex Over Time for Different Valsartan Doses');
% legend('0 mg','40 mg','80 mg','160 mg','Threshold');
% grid on;
% 
% figure('Name','Ang II-Receptor Complex for Different Doses','NumberTitle','off');
% 
% hold on;
% for i = 1:length(dose_levels)
%     plot(time_storage{i}, y_storage{i}(:,5), 'LineWidth', 2);  % Column 5: Ang II-Receptor Complex
% end
% 
% hold off;
% 
% xlabel('Time (hours)');
% ylabel('Concentration');
% title('Ang II-Receptor Complex Over Time for Different Valsartan Doses');
% legend('0 mg','40 mg','80 mg','160 mg','Threshold');
% grid on;
% 
% 
% figure('Name','All State Variables Over Time','NumberTitle','off');
% 
% subplot(2,2,1);  % 3x4 grid of subplots
% plot(time, y(:,1), 'LineWidth', 2);
% title(labels{1});
% xlabel('Time (hours)');
% ylabel('Concentration or Amount');
% grid on;
% 
% subplot(2,2,2);  % 3x4 grid of subplots
% plot(time, y(:,3), 'LineWidth', 2);
% title(labels{3});
% xlabel('Time (hours)');
% ylabel('Concentration or Amount');
% grid on;
% 
% subplot(2,2,3);  % 3x4 grid of subplots
% plot(time, y(:,4), 'LineWidth', 2);
% title(labels{4});
% xlabel('Time (hours)');
% ylabel('Concentration or Amount');
% grid on;
% 
% subplot(2,2,4);  % 3x4 grid of subplots
% plot(time, y(:,5), time, 8.34*10^(-5) * ones(size(time)),'LineWidth', 2);
% title(labels{5});
% xlabel('Time (hours)');
% ylabel('Concentration or Amount');
% grid on;
% 
% 
% sgtitle('Simulation Output: All State Variables');
