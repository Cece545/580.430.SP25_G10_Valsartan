function [out1,out2] = sim0_v2(p) %out1 is the time series, out2 is the 
% Define initial parameters
p.C0_2; %Resting concentration of total Receptors
p.C0_3; %Resting concentration of total Angiotensin II
p.doses; %List of tuples (time, dose) in h and mg
p.Vd = 17; %Volume of distribution L

% Define dosing times (hours) and amounts (mg)
doses = p.doses;

% Define time span for simulation
t_start = 0;
t_end   = max(p.doses(:,1)) + 24;  % simulate long enough after the last dose
time    = linspace(t_start, t_end, 1000);  % More refined time points
t_eval = linspace(t_start, t_end, 1000);

% Initial conditions:
y0 = [0, p.C0_2, p.C0_3, 0, 0, 0, 0, 0, 0, 0];  % Note: initial dose accounted for in doses, not here. Assuming no complexation at all.

% Initialize results storage
results_time = [];
results_concentration = [];
for i = 1:size(doses,1)
    dose_time  = doses(i,1);
    dose_mg    = doses(i,2);  % in mg
    dose_amt   = doses(i,2)*1e6/p.Vd/435.519 * 0.25;  % µM

    y0(10) = y0(10) + dose_amt;
% Set feedback ON or OFF based on time
    if dose_time < 336  % 2 weeks = 336 hours
        p.k_feedback = 0.000000984 / 7 * dose_mg * p.a;
    elseif dose_time < 672  % 4 weeks = 672 hours
        p.k_feedback = 0.000000984 / 7/3 * dose_mg * p.a;
    else
        p.k_feedback = 0;  % turn off feedback after 4 weeks
    end    % Define time span for this dose onward

    if i < size(doses,1)
        t_start_i = doses(i,1);
        t_end_i = doses(i+1,1);
    else
        t_start_i = doses(i,1);
        t_end_i = t_end;
    end

    % Skip if no time interval
    if t_start_i == t_end_i
        continue;
    end

    t_span = [t_start_i, t_end_i];
    options = odeset('RelTol',1e-6,'AbsTol',1e-9);
    [t_sol, y_sol] = ode15s(@(t,y) eqns_v2(t, y, p), t_span, y0, options);

    results_time = [results_time; t_sol];
    results_concentration = [results_concentration; y_sol];

    y0 = y_sol(end, :);
end

% Sort results to ensure correct time sequence
[results_time, idx] = sort(results_time);
results_concentration = results_concentration(idx, :);

out1 = results_time;
out2 = results_concentration;
end
