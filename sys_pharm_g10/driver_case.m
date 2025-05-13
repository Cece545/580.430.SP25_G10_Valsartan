close all; clear;

%% User-defined run case
RunCase = 7;  % Set to 1 or 2

%% Set global parameters
dose_interval = 24;  % hours
C0_2 = 0.08585;       % nM
C0_3 = 0.0239;        % nM
Vd = 17;              % L
p.a = 1;
switch RunCase
    case 1  % Regular on-time dosing
        
        p.C0_2 = C0_2;
        p.C0_3 = C0_3;
        p.Vd = Vd;
        dose_amount = 80;       % mg
        p.doses = [0, dose_amount];

        [t, y] = sim0_v2(p);
        y_receptor = (y(:, 1)+y(:, 6))*435.519/1000000;

        figure;
        plot(t, y_receptor, 'LineWidth', 2);
        title('Free Valsartan Concentration (On-Time Dosing)');
        xlabel('Time (h)');
        ylabel('Concentration (mg/L)');
        grid on;
    case 2  % Missed and re-taken dose at different times (generalized)

        % Base dosing schedule
        num_regular_doses = 4;   % First dose + scheduled second dose
        dose_interval = 24;     % Interval between doses (hours)
        dose_amount = 80;       % mg

        % Recovery time scenarios (relative to 2nd scheduled dose)
        recovery_offsets = dose_interval * (0:0.2:1);  % m/5 to m
        delays = (num_regular_doses - 1) * dose_interval + recovery_offsets;  % scheduled time + offset
        labels = ["Perfect", "m/5", "2m/5", "3m/5", "4m/5", "m"];
        cmap = lines(length(delays));  % color map

        % Initialize metrics table
        metrics = table('Size', [length(delays), 4], ...
                        'VariableTypes', ["string", "double", "double", "double"], ...
                        'VariableNames', ["Scenario", "Cmax", "Cmin", "AUC"]);

        figure('Name', 'Free Valsartan Concentration - All Scenarios');
        hold on;

        for i = 1:length(delays)
            % Build dose schedule
            recovery_time = delays(i);  % Time of missed/re-taken dose
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;

            % Regular doses up to (but not including) the missed one
            regular_times = (0:(num_regular_doses - 2)) * dose_interval;
            p.doses = [regular_times', repmat(dose_amount, length(regular_times), 1)];

            % Add the recovered dose
            p.doses = [p.doses; recovery_time, dose_amount];

            % Run simulation
            [t, y] = sim0_v2(p);
            D = y(:, 1);  % Free valsartan

            % Calculate metrics
            metrics.Scenario(i) = labels(i);
   
            % Trough is minimum in 24h before recovery (i.e. window before re-dose)
            trough_window = (recovery_time - dose_interval <= t) & (t <= recovery_time);
            auc_window = (t >= recovery_time) & (t <= recovery_time + 24);
            metrics.Cmax(i) = max(D(auc_window));
            metrics.Cmin(i) = min(D(auc_window));
            metrics.AUC(i) = trapz(t(auc_window), D(auc_window));

            % Plot concentration
            plot(t, D, 'LineWidth', 2, 'Color', cmap(i,:));
        end

        xlabel('Time (h)');
        ylabel('Free Valsartan (µM)');
        title('Free Valsartan Concentration vs. Time');
        legend(labels, 'Location', 'best');
        grid on;
        hold off;

        % 📊 Bar plots for PK metrics
        figure;
        subplot(1,3,1); bar(metrics.Cmax); title('C_{max}'); xticklabels(labels); ylabel('nM');
        subplot(1,3,2); bar(metrics.Cmin); title('C_{min}'); xticklabels(labels); ylabel('nM');
        subplot(1,3,3); bar(metrics.AUC);  title('AUC');     xticklabels(labels); ylabel('nM·h');
        sgtitle('Missed Dose Recovery - PK Metrics');

    case 3  % Missed and re-taken dose at different times

        delays = 24+dose_interval * (0:0.2:1);  % 0h to 24h (m)
        labels = ["Perfect", "m/5", "2m/5", "3m/5", "4m/5", "m"];
        cmap = lines(length(delays));  % color map
        metrics = table('Size',[length(delays),4], ...
                        'VariableTypes',["string","double","double","double"], ...
                        'VariableNames',["Scenario","Cmax","Cmin","AUC"]);
        figure('Name','Free Valsartan Concentration - All Scenarios');
        hold on;
        for i = 1:length(delays)
            recovery_time = delays(i);  % second dose time
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;
            p.doses = [0, dose_amount; recovery_time, dose_amount];

            [t, y] = sim0_v2(p);
            D = y(:,1);  % free valsartan

            metrics.Scenario(i) = labels(i);
            metrics.Cmax(i) = max(D);
            metrics.Cmin(i) = min(D(t >= recovery_time-24 & t <= recovery_time));  % trough before redosing
            metrics.AUC(i) = trapz(t, D);
            plot(t, D, 'LineWidth', 2, 'Color', cmap(i,:));

        end

        disp(metrics);

          % Finalize plot
        xlabel('Time (h)');
        ylabel('Free Valsartan (µM)');
        title('Free Valsartan Concentration vs. Time');
        legend(labels, 'Location', 'best');
        grid on;
        hold off;

        % 📊 Bar plots for PK metrics
        figure;
        subplot(1,3,1); bar(metrics.Cmax); title('C_{max}'); xticklabels(labels); ylabel('µM');
        subplot(1,3,2); bar(metrics.Cmin); title('C_{min}'); xticklabels(labels); ylabel('µM');
        subplot(1,3,3); bar(metrics.AUC);  title('AUC');     xticklabels(labels); ylabel('µM·h');
        sgtitle('Missed Dose Recovery - PK Metrics');
    case 4  % Missed second dose re-taken late + scheduled third dose at 48h

        dose_amount = 80;       % mg
        m = 24;                 % Standard dosing interval (hours)
        recovery_offsets = m * (0:5)/5;  % m/5 to m after 24h
        recovery_times = 24 + recovery_offsets;  % Time of re-taken second dose
        labels = ["0","m/5", "2m/5", "3m/5", "4m/5", "m"];
        cmap = lines(length(recovery_offsets));

        metrics = table('Size', [length(recovery_offsets), 4], ...
                        'VariableTypes', ["string", "double", "double", "double"], ...
                        'VariableNames', ["Scenario", "Cmax", "Cmin", "AUC"]);

        figure('Name','Free Valsartan Concentration - Missed 2nd Dose Recovered Late + 3rd Dose');
        hold on;

        for i = 1:length(recovery_times)
            recovery_time = recovery_times(i);  % Re-taken second dose
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;

            % Three doses: 0h, delayed second, and scheduled third at 48h
            p.doses = [0, dose_amount;
                       recovery_time, dose_amount;
                       48, dose_amount;
                       72, dose_amount;
                       96, dose_amount];

            % Run simulation
            [t, y] = sim0_v2(p);
            D = y(:,1);  % free valsartan

            % Calculate PK metrics
            metrics.Scenario(i) = labels(i);
            metrics.Cmax(i) = max(D);

            % Trough before 3rd dose (48h)
            % Trough before 3rd dose (48h), only if valid time points exis
            % Trough before 3rd dose (48h), only if valid time points exist
             % Trough before 3rd dose (48h), only if valid time points exist
            trough_window = (t >= recovery_time) & (t < 48);
            if any(trough_window)
                metrics.Cmin(i) = min(D(trough_window));
            else
                metrics.Cmin(i) = NaN;  % fallback if time window is empty
            end
            % AUC: 24h window after the recovery dose
            auc_window = (t >= recovery_time) & (t <= recovery_time + 24);
            metrics.AUC(i) = trapz(t(auc_window), D(auc_window));

            % Plot profile
            plot(t, D, 'LineWidth', 2, 'Color', cmap(i,:));
        end
        xlabel('Time (h)');
        ylabel('Free Valsartan (nM)');
        title('Free Valsartan Concentration Over Time');
        legend(labels, 'Location', 'best');
        grid on;
        hold off;

        % New Figure: Valsartan–ATR1 Complex Only
        figure('Name','Valsartan–ATR1 Complex Concentration');
        hold on;
        
        for i = 1:length(recovery_times)
            recovery_time = recovery_times(i);
        
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;
        % Three doses: 0h, delayed second, and scheduled third at 48h
            p.doses = [0, dose_amount;
                       recovery_time, dose_amount;
                       48, dose_amount;
                       72, dose_amount;
                       96, dose_amount];
            [t, y] = sim0_v2(p);
            DR = y(:, 5);  % Valsartan–ATR1 complex
        
            plot(t, DR, 'LineWidth', 2, 'Color', cmap(i,:));
        end
        
        xlabel('Time (h)');
        ylabel('AngII–ATR1 Complex (nM)');
        title('AngII–ATR1 Complex Concentration Over Time');
        legend(labels, 'Location', 'best');
        grid on;
        hold off;


        % Plot metrics
        figure;
        subplot(1,3,1); bar(metrics.Cmax); title('C_{max}'); xticklabels(labels); ylabel('µM');
        subplot(1,3,2); bar(metrics.Cmin); title('C_{min}'); xticklabels(labels); ylabel('µM');
        subplot(1,3,3); bar(metrics.AUC);  title('AUC (24h after Recovery)'); xticklabels(labels); ylabel('µM·h');
        sgtitle('Impact of Late Re-Taken Dose on PK');
      case 5  % 🟢 Simulate different normal doses: 0mg, 40mg, 80mg, 160mg

        dose_levels = [0, 40, 80, 160];  % mg
        dose_interval = 24;             % every 24h
        num_doses = 10;                  % number of doses (e.g. 0h, 24h, 48h)
        labels = ["0 mg", "40 mg", "80 mg", "160 mg"];
        % Generate cold shades of blue (light to dark)
        num_colors = length(dose_levels);
        cmap = zeros(num_colors, 3);
        cmap = [
            0.8, 0.9, 1.0;   % light icy blue (0 mg)
            0.4, 0.7, 1.0;   % medium blue (40 mg)
            0.2, 0.5, 0.9;   % strong mid-blue (80 mg)
            0.0, 0.2, 0.6    % deep navy blue (160 mg)
        ];        
        figure('Name','Free Valsartan for Different Doses');
        hold on;

        figure('Name','Valsartan–ATR1 Complex for Different Doses');
        hold on;
        figure('Name','Angiotensin II–ATR1 Complex for Different Doses');
        hold on;

        for i = 1:length(dose_levels)
            dose = dose_levels(i);

            % Build dose schedule
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;

            times = (0:num_doses-1) * dose_interval;
            p.doses = [times', repmat(dose, num_doses, 1)];

            % Run simulation
            [t, y] = sim0_v2(p);
            D = y(:,1);   % Free valsartan
            DR = y(:,4);  % Valsartan–ATR1 complex
            AR = y(:,5);  % Angeotensin–ATR1 complex

            % Plot Free Valsartan
            figure(1);
            plot(t, D, 'LineWidth', 2, 'Color', cmap(i,:));

            % Plot Complex
            figure(2);
            plot(t, DR, 'LineWidth', 2, 'Color', cmap(i,:));
            % Plot Complex
            figure(3);
            plot(t, AR, 'LineWidth', 2, 'Color', cmap(i,:));
        end

        % Finalize Free Drug Plot
        figure(1);
        xlabel('Time (h)');
        ylabel('Free Valsartan (nM)');
        title('Free Valsartan for Different Doses');
        legend(labels, 'Location', 'northeast');
        grid on; hold off;

        % Finalize Complex Plot
        figure(2);
        xlabel('Time (h)');
        ylabel('Valsartan–ATR1 Complex (nM)');
        title('Valsartan–ATR1 Complex for Different Doses');
        legend(labels, 'Location', 'northeast');
        grid on; hold off;
        % 🧬 Plot Angiotensin II concentration (y(:,3))
        figure('Name','Angiotensin II for Different Doses');
        hold on;

        % Finalize Complex Plot
        figure(3);
        xlabel('Time (h)');
        ylabel('Angiotensin II–ATR1 Complex (nM)');
        title('Angiotensin II–ATR1 Complex for Different Doses');
        legend(labels, 'Location', 'northeast');
        grid on; hold off;
        % 🧬 Plot Angiotensin II concentration (y(:,3))
        figure('Name','Angiotensin II for Different Doses');
        hold on;
        for i = 1:length(dose_levels)
            dose = dose_levels(i);
        
            % Reconstruct p and run sim again to access y(:,3)
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;
        
            times = (0:num_doses-1) * dose_interval;
            p.doses = [times', repmat(dose, num_doses, 1)];
        
            [t, y] = sim0_v2(p);
            AII = y(:,3);  % Angiotensin II
        
            plot(t, AII, 'LineWidth', 2.5, 'Color', cmap(i,:));
        end
        
        xlabel('Time (h)');
        ylabel('Angiotensin II (nM)');
        title('Angiotensin II Concentration for Different Doses');
        legend(labels, 'Location', 'northeast');
        grid on;
        hold off;
% 🧪 Compare with and without feedback after missing a dose
    case 6
    
        dose_interval = 24;         % hours between doses
        dose_amount = 80;           % mg per dose
        num_total_doses = 10;       % total scheduled doses
        miss_dose_idx = [5];     % doses to skip;       
        
        % Create time schedule and mask out missed dose
        dose_times = (0:num_total_doses-1) * dose_interval;
        dose_mask = true(size(dose_times));
        dose_mask(miss_dose_idx) = false;
        
        % Common dose schedule (excluding missed one)
        active_times = dose_times(dose_mask);
        doses_base = [active_times', repmat(dose_amount, numel(active_times), 1)];
        
        % Preallocate color and labels
        cmap = [0.3 0.6 1.0; 0.1 0.1 0.5; 0.6 0.6 0.6];  % light blue, dark blue, gray
        labels = ["With Feedback", "Without Feedback", "No Drug"];
        
        figure('Name','Angiotensin II–Receptor Complex With vs Without Feedback');
        hold on;
        
        for k = 1:3
            p.C0_2 = C0_2;
            p.C0_3 = C0_3;
            p.Vd = Vd;
            
            if k == 1
                % With feedback
                p.a = 1;
                p.doses = doses_base;
                line_style = '--';
    
            elseif k == 2
                % No feedback
                p.a = 0;
                p.doses = doses_base;
                line_style = '-';
    
            else
                % No drug
                p.a = 1;  % feedback doesn't matter — no drug given
                p.doses = [dose_times', zeros(num_total_doses, 1)];
                line_style = ':';
            end
    
            [t, y] = sim0_v2(p);
            AR = y(:,5);  % Angiotensin II–Receptor Complex
            plot(t, AR, 'LineWidth', 2.5, 'Color', cmap(k,:), 'LineStyle', line_style);
        end
    
        xlabel('Time (h)');
        ylabel('Ang II–Receptor Complex (µM)');
        title('Ang II–Receptor Complex with Doses #5 Missed');
        legend(labels, 'Location', 'northeast');
        grid on;
        hold off;
    case 7
        dose_amount = 80;   % mg
        m = 24;             % standard dosing interval (hours)
        missed_doses = [2, 3, 4, 5, 6];  % indices of the dose to miss and delay
        recovery_offsets = m * (0:5)/5;
        recovery_labels = ["0", "m/5", "2m/5", "3m/5", "4m/5", "m"];
        cmap = lines(length(recovery_offsets));
        
        for miss_idx = 1:length(missed_doses)
            missed_dose_number = missed_doses(miss_idx);
            
            metrics = table('Size', [length(recovery_offsets), 5], ...
                            'VariableTypes', ["string", "double", "double", "double", "double"], ...
                            'VariableNames', ["Scenario", "MissedDose", "Cmax", "Cmin", "AUC"]);
        
            figure('Name', sprintf('Free Valsartan - Missed Dose #%d', missed_dose_number));
            hold on;
        
            for i = 1:length(recovery_offsets)
                recovery_time = m * (missed_dose_number - 1) + recovery_offsets(i);  % Delayed dose time
                
                % Build doses, omitting the missed dose
                dose_schedule = [];
                for d = 0:7
                    t_dose = m * d;
                    if d + 1 ~= missed_dose_number  % skip missed dose
                        dose_schedule = [dose_schedule; t_dose, dose_amount];
                    end
                end
                % Add delayed missed dose
                dose_schedule = [dose_schedule; recovery_time, dose_amount];
                % Sort doses in chronological order
                dose_schedule = sortrows(dose_schedule, 1);
        
                % Assign to parameter struct
                p.doses = dose_schedule;
                p.C0_2 = C0_2;
                p.C0_3 = C0_3;
                p.Vd = Vd;
        
                % Simulate
                [t, y] = sim0_v2(p);
                D = y(:, 1);       % free valsartan
                AR = y(:, 5);      % Ang II–ATR1 complex
        
                % Calculate metrics
                metrics.Scenario(i) = recovery_labels(i);
                metrics.MissedDose(i) = missed_dose_number;
                metrics.Cmax(i) = max(D);
        
                % Cmin before next scheduled dose after recovery
                next_dose_time = min(dose_schedule(dose_schedule(:,1) > recovery_time, 1), [], 'omitnan');
                if ~isempty(next_dose_time)
                    trough_window = (t >= recovery_time) & (t < next_dose_time);
                else
                    trough_window = (t >= recovery_time);
                end
                metrics.Cmin(i) = min(D(trough_window));
        
                % AUC: 24h window after recovery
                auc_window = (t >= recovery_time) & (t <= recovery_time + 24);
                metrics.AUC(i) = trapz(t(auc_window), D(auc_window));
        
                % Plot Free Valsartan
                plot(t, D, 'LineWidth', 2, 'Color', cmap(i, :));

                % Save time series for Free Valsartan and AngII–ATR1 Complex
                time_series_table = table(t, D, AR, ...
                    'VariableNames', {'Time_hr', 'Free_Valsartan_nM', 'AngII_ATR1_Complex_nM'});

                safe_label = strrep(recovery_labels(i), '/', '_');  % Replace '/' with '_'
                filename_ts = sprintf('missed_dose_%d_recovery_%s_timeseries.csv', ...
                missed_dose_number, safe_label);

                writetable(time_series_table, filename_ts);
            end
        
            xlabel('Time (h)');
            ylabel('Free Valsartan (nM)');
            title(sprintf('Free Valsartan - Missed Dose #%d Recovered Late', missed_dose_number));
            legend(recovery_labels, 'Location', 'best');
            grid on;
            hold off;
        
            % Save metrics table for this scenario
            filename_csv = sprintf('missed_dose_%d_metrics.csv', missed_dose_number);
            writetable(metrics, filename_csv);
        
            % Optional: save .mat too
            filename_mat = sprintf('missed_dose_%d_metrics.mat', missed_dose_number);
            save(filename_mat, 'metrics');
        end

    case 8  % Extract AngII–ATR1 Complex values at select times for different doses
        dose_levels = [10, 40, 80, 160];  % mg
        timepoints = [0, 2, 4, 6, 336, 338, 340, 342, 672, 674, 676, 678];  % hours
        Vd = 17;
        p.Vd = Vd;
        p.C0_2 = C0_2;
        p.C0_3 = C0_3;

        % Table to collect results
        all_results = table();

        for i = 1:length(dose_levels)
            dose = dose_levels(i);
            % Build dose schedule: simulate over 5 weeks
            dose_interval = 24;
            num_doses = 30;  % enough to cover ~700 hours
            times = (0:num_doses-1) * dose_interval;
            p.doses = [times', repmat(dose, num_doses, 1)];

            % Run simulation
            [t, y] = sim0_v2(p);
            AR = y(:,5);  % Ang II–ATR1 Complex

            % Remove duplicate time entries before interpolation
            [t_unique, unique_idx] = unique(t);
            AR_unique = AR(unique_idx);

            % Interpolate AR values at specified timepoints
            try
                AR_interp = interp1(t_unique, AR_unique, timepoints, 'linear', 'extrap');
            catch
                warning('Interpolation failed for dose %d. Filling with NaNs.', dose);
                AR_interp = nan(size(timepoints));
            end

            % Ensure vector lengths match
            if length(AR_interp) ~= length(timepoints)
                warning('Interpolated vector length mismatch at dose %d', dose);
                AR_interp = nan(size(timepoints));
            end

            % Store in table
            dose_col = repmat(dose, length(timepoints), 1);
            time_col = timepoints(:);  % force column vector
            result_table = table(dose_col, time_col, AR_interp(:), ...
                            'VariableNames', {'Dose_mg', 'Time_hr', 'AngII_ATR1_Complex_nM'});

            all_results = [all_results; result_table];
        end

        % Display and save
        disp(all_results);
        writetable(all_results, 'AngII_ATR1_Complex_Timepoints.csv');


end

