clc; clear; close all;

%% Step 1: Load and Sample 500 Random Values from Column 2
% Load all data from column 2
traditional_all = readmatrix('normalvolt_data.csv', 'Range', 'B:B');
selfhealing_all = readmatrix('voltage_varing_data.csv', 'Range', 'B:B');

% Remove NaN values
traditional_all = traditional_all(~isnan(traditional_all));
selfhealing_all = selfhealing_all(~isnan(selfhealing_all));

% Check if files have enough data
min_samples = 500;
if length(traditional_all) < min_samples || length(selfhealing_all) < min_samples
    error('Error: One or both files contain less than 500 data points.');
end

% Randomly sample 500 points (without replacement)
rng(250); % For reproducibility
traditional_data = datasample(traditional_all, min_samples, 'Replace', false);
selfhealing_data = datasample(selfhealing_all, min_samples, 'Replace', false);

%% Step 2: Compute Descriptive Statistics (unchanged)
mu_traditional = mean(traditional_data);
mu_selfhealing = mean(selfhealing_data);

std_traditional = std(traditional_data);
std_selfhealing = std(selfhealing_data);

fprintf('Traditional System: Mean = %.2f v, Std Dev = %.2f v (n=%d)\n', ...
        mu_traditional, std_traditional, length(traditional_data));
fprintf('Self-Healing System: Mean = %.2f v, Std Dev = %.2f v (n=%d)\n\n', ...
        mu_selfhealing, std_selfhealing, length(selfhealing_data));

%% [Rest of the code remains identical from Step 3 onward]
%% Step 3: Perform Two-Sample t-Test
alpha = 0.5;
[h, p, ci, stats] = ttest2(selfhealing_data, traditional_data, 'Alpha', alpha, 'Tail', 'both');

fprintf('--- Hypothesis Test Results ---\n');
fprintf('t-statistic = %.4f\n', stats.tstat);
fprintf('Degrees of freedom = %.2f\n', stats.df);
fprintf('p-value = %.6f\n', p);

if h == 1
    fprintf('\nConclusion: REJECT H0 (Self-healing is significantly faster)\n');
else
    fprintf('\nConclusion: FAIL TO REJECT H0 (No significant difference)\n');
end

%% Step 4: Effect Size
pooled_std = sqrt((std_traditional^2 + std_selfhealing^2)/2);
cohen_d = (mu_traditional - mu_selfhealing)/pooled_std;
fprintf('Effect Size (Cohen''s d) = %.2f\n', cohen_d);

%% Step 5: Plotting
figure;
boxplot([traditional_data, selfhealing_data], 'Labels', {'Traditional', 'Self-Healing'});
ylabel('Recovery voltage (v)');
title(sprintf('Comparison of Recovery Times\n500 Random Samples per System'));
hold on;
plot(1, mu_traditional, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
plot(2, mu_selfhealing, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;
grid on;

%% Step 6: Save Results
results = struct(...
    'Traditional_Mean', mu_traditional, ...
    'SelfHealing_Mean', mu_selfhealing, ...
    'Traditional_Std', std_traditional, ...
    'SelfHealing_Std', std_selfhealing, ...
    't_statistic', stats.tstat, ...
    'df', stats.df, ...
    'p_value', p, ...
    'Cohen_d', cohen_d, ...
    'SampleSize', min_samples);
save('hypothesis_test_results.mat', 'results');
