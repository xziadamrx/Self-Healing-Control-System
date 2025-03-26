%% Clear environment and setup
clc; clearvars;

%% 1) Simulation Setup
T = 1000;            % Total simulation steps
dt = 0.001;          % Time step
time = 0:dt:(T*dt);  % Time vector
nSteps = length(time);

% Normal motor values
normal_params.voltage     = 400 * sin(2*pi*50*time);
normal_params.current     = 10  * sin(2*pi*50*time);
normal_params.temperature = 75;
normal_params.vibration   = 2;

% Standard deviations for variations
std_dev.voltage     = 20;
std_dev.current     = 2;
std_dev.temperature = 5;
std_dev.vibration   = 1;

% Prior probabilities
P_normal = 0.9;
P_fault  = 0.1;

% Monte Carlo parameters
num_simulations = 1000;
fault_prob_results = zeros(num_simulations, nSteps);

for sim = 1:num_simulations
    % Initialize arrays
    voltage_arr       = normal_params.voltage;
    current_arr       = normal_params.current;
    temperature_arr   = normal_params.temperature + randn(1, nSteps)*std_dev.temperature;
    vibration_arr     = normal_params.vibration + randn(1, nSteps)*std_dev.vibration;
    fault_prob_arr    = 0.1 * ones(1, nSteps);
    final_voltage_arr = normal_params.voltage;
    final_current_arr = normal_params.current;

    %% 2) Introduce Faults & Apply Adjustments
    for k = floor(nSteps/2):nSteps
        % Apply a random fault scale factor
        fault_scale = 1 + 0.3 * randn();
        voltage_arr(k) = normal_params.voltage(k) * fault_scale;
        current_arr(k) = normal_params.current(k) * fault_scale;
        
        % Initialize correction factor
        correction_factor = 1;
        
        % Apply correction if voltage or current exceed thresholds
        if abs(voltage_arr(k)) > abs(normal_params.voltage(k)) * 1.05 || ...
           abs(current_arr(k)) > abs(normal_params.current(k)) * 1.1
            correction_factor = correction_factor * 0.8;
        end
        
        if temperature_arr(k) > 80
        fprintf('Warning: High temperature at %.2f s!\n', time(k));
        correction_factor = correction_factor * 0.8; 
    end

    % Check for high vibration fault and adjust the correction factor
    if vibration_arr(k) > 3
        fprintf('Warning: High vibration at %.2f s!\n', time(k));
        correction_factor = correction_factor * 0.8; 
    end

        
        % Prevent excessive correction
        correction_factor = max(correction_factor, 0.5);
        
        % Apply final correction
        final_voltage_arr(k) = normal_params.voltage(k) * correction_factor;
        final_current_arr(k) = normal_params.current(k) * correction_factor;
    end

    %% 3) Compute Overall Fault Probability using Bayesian Inference
    for k = floor(nSteps/2):nSteps
        measured.voltage     = final_voltage_arr(k);
        measured.current     = final_current_arr(k);
        measured.temperature = temperature_arr(k);
        measured.vibration   = vibration_arr(k);
        
        fault_prob_arr(k) = compute_fault_probability(measured, normal_params, std_dev, P_normal, P_fault, k);
    end
    
    % Store results for Monte Carlo analysis
    fault_prob_results(sim, :) = fault_prob_arr;
end

% Compute mean fault probability across simulations
mean_fault_prob = mean(fault_prob_results, 1);

% Plot Monte Carlo fault probability
figure;
plot(time, mean_fault_prob, 'r', 'LineWidth', 1.5);
title('Monte Carlo Estimated Fault Probability');
xlabel('Time (s)');
ylabel('Probability');
grid on;
%% 4) Subplots for Data and Fault Probability
figure('Name','Motor Simulation Results','NumberTitle','off');

% Subplot 1: Normal Voltage vs Time (Blue)
subplot(4,2,1)
plot(time, normal_params.voltage, 'b');
title('Normal Voltage');
xlabel('Time (s)');
ylabel('Voltage (V)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(normal_params.voltage) - min(normal_params.voltage));
ylim([min(normal_params.voltage)-margin, max(normal_params.voltage)+margin]);

% Subplot 2: Normal Current vs Time (Red)
subplot(4,2,2)
plot(time, normal_params.current, 'r');
title('Normal Current');
xlabel('Time (s)');
ylabel('Current (A)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(normal_params.current) - min(normal_params.current));
ylim([min(normal_params.current)-margin, max(normal_params.current)+margin]);

% Subplot 3: Voltage Array vs Time (Green)
subplot(4,2,3)
plot(time, voltage_arr, 'g');
title('Voltage Array');
xlabel('Time (s)');
ylabel('Voltage (V)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(voltage_arr) - min(voltage_arr));
ylim([min(voltage_arr)-margin, max(voltage_arr)+margin]);

% Subplot 4: Current Array vs Time (Magenta)
subplot(4,2,4)
plot(time, current_arr, 'm');
title('Current Array');
xlabel('Time (s)');
ylabel('Current (A)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(current_arr) - min(current_arr));
ylim([min(current_arr)-margin, max(current_arr)+margin]);

% Subplot 5: Final Voltage Array vs Time (Cyan)
subplot(4,2,5)
plot(time, final_voltage_arr, 'c');
title('Final Voltage Array');
xlabel('Time (s)');
ylabel('Voltage (V)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(final_voltage_arr) - min(final_voltage_arr));
ylim([min(final_voltage_arr)-margin, max(final_voltage_arr)+margin]);

% Subplot 6: Final Current Array vs Time (Black)
subplot(4,2,6)
plot(time, final_current_arr, 'k');
title('Final Current Array');
xlabel('Time (s)');
ylabel('Current (A)');
xlim([time(1) time(end)]);
margin = 0.2 * (max(final_current_arr) - min(final_current_arr));
ylim([min(final_current_arr)-margin, max(final_current_arr)+margin]);

% Subplot 7 (spanning two columns): Fault Probability vs Time (Dark Orange)
subplot(4,2,[7 8])
plot(time, fault_prob_arr, 'Color',[0.85 0.33 0.10],'LineWidth',1.5);
title('Fault Probability');
xlabel('Time (s)');
ylabel('Probability');
xlim([time(1) time(end)]);
ylim([0 1]);
grid on;
%% Bayesian Fault Probability Function
function prob = compute_fault_probability(measured, normal_params, std_dev, P_normal, P_fault, time_idx)
    deviation = abs([ ...
        measured.voltage      - normal_params.voltage(time_idx), ...
        measured.current      - normal_params.current(time_idx), ...
        measured.temperature  - normal_params.temperature,       ...
        measured.vibration    - normal_params.vibration          ...
    ]);
    
    thresholds = [std_dev.voltage, std_dev.current, std_dev.temperature, std_dev.vibration] ;
    
    likelihood_normal_total = 1;
    for i = 1:4
        param_likelihood = (1/sqrt(2*pi)) * exp(-0.5 * (deviation(i)/thresholds(i))^2);
        likelihood_normal_total = likelihood_normal_total * param_likelihood;
    end
    
    likelihood_fault_total = 1 - likelihood_normal_total;
    
    P_normal_post = (likelihood_normal_total * P_normal) / ...
                    (likelihood_normal_total * P_normal + likelihood_fault_total * P_fault);
    
    P_fault_post = 1 - P_normal_post;
    
    prob = min(1, max(0, P_fault_post));
end
