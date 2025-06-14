function results = exp1_pitch_modification()
% EXP1_PITCH_MODIFICATION - Experiment 1: Pitch Modification with Periodic Signals
%
% This experiment generates a periodic signal and applies time-scale 
% modification to transpose it to different frequencies.
%
% Output:
%   results - Structure containing experiment results

fprintf('=== Experiment 1: Pitch Modification with Periodic Signals ===\n\n');

% Get configuration
config = project_config('exp1');

% Step 1: Generate original signal (C4 + G5)
fprintf('Step 1: Generating original signal...\n');
[x_original, t] = generate_periodic_signal(config.duration);

% Step 2: Apply time-scale modification to transpose to target frequencies
fprintf('\nStep 2: Applying time-scale modification...\n');

% Calculate stretch factors to achieve target frequencies
% To transpose from C4 (261 Hz) to G4 (392 Hz)
alpha_g4 = config.target_g4 / config.freq1;  % 392/261 ≈ 1.5
% To transpose from C4 (261 Hz) to F3 (175 Hz)  
alpha_f3 = config.target_f3 / config.freq1;  % 175/261 ≈ 0.67

% Apply OLA
y_ola_g4 = ola(x_original, alpha_g4);
y_ola_f3 = ola(x_original, alpha_f3);

% Apply WSOLA  
y_wsola_g4 = wsolaTSM(x_original, alpha_g4);
y_wsola_f3 = wsolaTSM(x_original, alpha_f3);

fprintf('Applied TSM: α_G4 = %.2f, α_F3 = %.2f\n', alpha_g4, alpha_f3);

% Step 3: Plot results
fprintf('\nStep 3: Plotting amplitude vs time...\n');

% Create output directory for plots
plots_dir = fullfile('outputs', 'experiment1', 'plots');
mkdir(plots_dir);

% Plot comparison of original vs OLA vs WSOLA for G4 transposition
signals_g4 = {x_original, y_ola_g4, y_wsola_g4};
titles_g4 = {'Original (C4+G5)', sprintf('OLA → G4 (α=%.2f)', alpha_g4), sprintf('WSOLA → G4 (α=%.2f)', alpha_g4)};
plot_amplitude_vs_time(signals_g4, t, titles_g4, [0, 0.02]);  % First 20ms
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_G4_comparison.png'));

% Plot comparison for F3 transposition
signals_f3 = {x_original, y_ola_f3, y_wsola_f3};
titles_f3 = {'Original (C4+G5)', sprintf('OLA → F3 (α=%.2f)', alpha_f3), sprintf('WSOLA → F3 (α=%.2f)', alpha_f3)};
plot_amplitude_vs_time(signals_f3, t, titles_f3, [0, 0.02]);  % First 20ms
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_F3_comparison.png'));

fprintf('Plots saved to: %s\n', plots_dir);

% Step 4: Save audio files
fprintf('\nStep 4: Saving audio files...\n');

% Create output directory for experiment 1
output_dir = fullfile('outputs', 'experiment1', 'audio');

% Fixed sampling rate
fs = 44100;

% Save original signal
audiowrite(fullfile(output_dir, 'original_C4_G5.wav'), x_original, fs);

% Save OLA results
audiowrite(fullfile(output_dir, sprintf('OLA_G4_alpha_%.2f.wav', alpha_g4)), y_ola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('OLA_F3_alpha_%.2f.wav', alpha_f3)), y_ola_f3, fs);

% Save WSOLA results
audiowrite(fullfile(output_dir, sprintf('WSOLA_G4_alpha_%.2f.wav', alpha_g4)), y_wsola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('WSOLA_F3_alpha_%.2f.wav', alpha_f3)), y_wsola_f3, fs);

fprintf('Audio files saved to: %s\n', output_dir);

% Step 5: Store results
results.config = config;
results.signals.original = x_original;
results.signals.ola_g4 = y_ola_g4;
results.signals.ola_f3 = y_ola_f3;
results.signals.wsola_g4 = y_wsola_g4;
results.signals.wsola_f3 = y_wsola_f3;
results.time = t;
results.stretch_factors.alpha_g4 = alpha_g4;
results.stretch_factors.alpha_f3 = alpha_f3;

fprintf('\nExperiment 1 completed successfully!\n');
fprintf('Results stored in output structure.\n');

end