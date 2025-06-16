function results = exp1_pitch_modification()
% EXP1_PITCH_MODIFICATION - Experiment 1: Pitch Modification with Periodic Signals
%
% This experiment generates a periodic signal and applies pitch shifting
% to transpose it to G4 (Dominant) and F3 (Sub-Dominant) frequencies.
%
% Output:
%   results - Structure containing experiment results

fprintf('=== Experiment 1: Pitch Modification with Periodic Signals ===\n\n');

% Fixed sampling rate
fs = 44100;

% Get configuration
config = project_config('exp1');

% Step 1: Generate original signal (C4 + G5)
fprintf('Step 1: Generating original signal...\n');
[x_original, t] = generate_periodic_signal(config.duration);

% Step 2: Apply pitch shifting to transpose to target frequencies
fprintf('\nStep 2: Applying pitch shifting...\n');

% Calculate pitch shift in cents
% From C4 (261 Hz) to G4 (392 Hz): ~700 cents
cents_g4 = 1200 * log2(config.target_g4 / config.freq1);
% From C4 (261 Hz) to F3 (175 Hz): ~-700 cents  
cents_f3 = 1200 * log2(config.target_f3 / config.freq1);

fprintf('Pitch shifts: G4 = %.0f cents, F3 = %.0f cents\n', cents_g4, cents_f3);

% Set up parameters for pitchShiftViaTSM
param_ola.fsAudio = fs;
param_ola.algTSM = @(x, alpha, ~) ola(x, alpha);

param_wsola.fsAudio = fs;
param_wsola.algTSM = @wsolaTSM;

% Apply pitch shifting
y_ola_g4 = pitchShiftViaTSM(x_original, cents_g4, param_ola);
y_ola_f3 = pitchShiftViaTSM(x_original, cents_f3, param_ola);
y_wsola_g4 = pitchShiftViaTSM(x_original, cents_g4, param_wsola);
y_wsola_f3 = pitchShiftViaTSM(x_original, cents_f3, param_wsola);

% Create plots directory for FFT plots
plots_dir = fullfile('outputs', 'experiment1', 'plots');

% Step 3: Plot amplitude vs time comparisons
fprintf('\nStep 3: Plotting amplitude vs time comparisons...\n');

% Plot comparison for G4 pitch shift
signals_g4 = {x_original, y_ola_g4, y_wsola_g4};
titles_g4 = {'Original (C4+G5)', sprintf('OLA → G4 (%.0f cents)', cents_g4), sprintf('WSOLA → G4 (%.0f cents)', cents_g4)};
plot_amplitude_vs_time(signals_g4, t, titles_g4, [0, 0.02]);
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_G4.png'));  % Add this line

% Plot comparison for F3 pitch shift
signals_f3 = {x_original, y_ola_f3, y_wsola_f3};
titles_f3 = {'Original (C4+G5)', sprintf('OLA → F3 (%.0f cents)', cents_f3), sprintf('WSOLA → F3 (%.0f cents)', cents_f3)};
plot_amplitude_vs_time(signals_f3, t, titles_f3, [0, 0.02]);
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_F3.png'));  % Add this line

% Step 4: FFT Analysis
fprintf('\nStep 4: Analyzing frequency content...\n');

% Analyze original signal
fprintf('\n--- ORIGINAL SIGNAL (C4 + G5) ---\n');
[freq_orig, mag_orig, peaks_orig] = analyze_fft(x_original, fs);
title('Original Signal FFT (C4+G5)');
saveas(gcf, fullfile(plots_dir, 'fft_original_signal.png'));

% Analyze OLA results
fprintf('\n--- OLA G4 PITCH SHIFT (%.0f cents) ---\n', cents_g4);
[freq_ola_g4, mag_ola_g4, peaks_ola_g4] = analyze_fft(y_ola_g4, fs);
title(sprintf('OLA G4 Pitch Shift FFT (%.0f cents)', cents_g4));
saveas(gcf, fullfile(plots_dir, 'fft_ola_g4.png'));

fprintf('\n--- OLA F3 PITCH SHIFT (%.0f cents) ---\n', cents_f3);
[freq_ola_f3, mag_ola_f3, peaks_ola_f3] = analyze_fft(y_ola_f3, fs);
title(sprintf('OLA F3 Pitch Shift FFT (%.0f cents)', cents_f3));
saveas(gcf, fullfile(plots_dir, 'fft_ola_f3.png'));

% Analyze WSOLA results
fprintf('\n--- WSOLA G4 PITCH SHIFT (%.0f cents) ---\n', cents_g4);
[freq_wsola_g4, mag_wsola_g4, peaks_wsola_g4] = analyze_fft(y_wsola_g4, fs);
title(sprintf('WSOLA G4 Pitch Shift FFT (%.0f cents)', cents_g4));
saveas(gcf, fullfile(plots_dir, 'fft_wsola_g4.png'));

fprintf('\n--- WSOLA F3 PITCH SHIFT (%.0f cents) ---\n', cents_f3);
[freq_wsola_f3, mag_wsola_f3, peaks_wsola_f3] = analyze_fft(y_wsola_f3, fs);
title(sprintf('WSOLA F3 Pitch Shift FFT (%.0f cents)', cents_f3));
saveas(gcf, fullfile(plots_dir, 'fft_wsola_f3.png'));

fprintf('FFT plots saved to: %s\n', plots_dir);

% Step 5: Save audio files
fprintf('\nStep 5: Saving audio files...\n');

output_dir = fullfile('outputs', 'experiment1', 'audio');

% Save all signals
audiowrite(fullfile(output_dir, 'original_C4_G5.wav'), x_original, fs);
audiowrite(fullfile(output_dir, sprintf('OLA_G4_%.0f_cents.wav', cents_g4)), y_ola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('OLA_F3_%.0f_cents.wav', cents_f3)), y_ola_f3, fs);
audiowrite(fullfile(output_dir, sprintf('WSOLA_G4_%.0f_cents.wav', cents_g4)), y_wsola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('WSOLA_F3_%.0f_cents.wav', cents_f3)), y_wsola_f3, fs);

fprintf('Audio files saved to: %s\n', output_dir);

% Step 6: Store results
results.config = config;
results.fs = fs;
results.signals.original = x_original;
results.signals.ola_g4 = y_ola_g4;
results.signals.ola_f3 = y_ola_f3;
results.signals.wsola_g4 = y_wsola_g4;
results.signals.wsola_f3 = y_wsola_f3;
results.time = t;
results.pitch_shifts.cents_g4 = cents_g4;
results.pitch_shifts.cents_f3 = cents_f3;
results.fft_analysis.original = peaks_orig;
results.fft_analysis.ola_g4 = peaks_ola_g4;
results.fft_analysis.ola_f3 = peaks_ola_f3;
results.fft_analysis.wsola_g4 = peaks_wsola_g4;
results.fft_analysis.wsola_f3 = peaks_wsola_f3;

fprintf('\nExperiment 1 completed successfully!\n');
fprintf('Results stored in output structure.\n');

end