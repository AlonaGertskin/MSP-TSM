function results = exp1_pitch_modification()
% EXP1_PITCH_MODIFICATION - Experiment 1: Pitch Modification with Periodic Signals
%
% This experiment systematically compares OLA and WSOLA algorithms for pitch 
% modification using a controlled periodic test signal. The experiment performs
% parameter optimization to find optimal settings for both algorithms, then
% applies pitch shifting to transpose the signal to musically relevant frequencies.
%
% Process:
%   1. Generates periodic test signal: x = cos(2π×261×t) + 0.4×cos(2π×783×t)
%      (C4 + G5 harmonic combination)
%   2. Performs parameter sweep to optimize OLA settings (synHop, winLen, winType)
%   3. Performs tolerance sweep to optimize WSOLA settings (using best OLA params)
%   4. Applies pitch shifting to transpose signal to:
%      - G4 (392 Hz, +704 cents) - Dominant
%      - F3 (175 Hz, -692 cents) - Sub-Dominant  
%   5. Generates comparative analysis plots (amplitude vs time, FFT spectra)
%   6. Saves optimized audio files and comprehensive parameter sweep documentation
%
% Output:
%   results - Structure containing:
%             .optimal_params - Best parameters found for OLA and WSOLA
%             .signals - All processed signals (original, OLA, WSOLA)
%             .fft_analysis - Frequency analysis results for all signals
%             .pitch_shifts - Cent values for target frequencies
%
% Files Generated:
%   - Audio: Original and pitch-shifted signals using optimal parameters
%   - Plots: Amplitude comparisons, FFT spectra, parameter sweep results
%   - Log: Complete parameter sweep results and optimization summary

fprintf('=== Experiment 1: Pitch Modification with Periodic Signals ===\n\n');

% Fixed sampling rate
fs = 44100;

% Get configuration
config = project_config('exp1');

% Step 1: Generate original signal (C4 + G5)
fprintf('Step 1: Generating original signal...\n');
[x_original, t] = generate_periodic_signal(config.duration);

% Create plots directory for FFT plots
plots_dir = fullfile('outputs', 'experiment1', 'plots');

% Step 2: Parameter sweep to find optimal settings
fprintf('\nStep 2: Parameter sweep to find optimal settings...\n');

% Calculate pitch shift in cents
% From C4 (261 Hz) to G4 (392 Hz): ~700 cents
cents_g4 = 1200 * log2(config.target_g4 / config.freq1);
% From C4 (261 Hz) to F3 (175 Hz): ~-700 cents  
cents_f3 = 1200 * log2(config.target_f3 / config.freq1);

% Create log file for parameter sweep results
log_file = fullfile('outputs', 'experiment1', 'parameter_sweep_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 1: Parameter Sweep Results ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));
fprintf(fid, 'Target frequencies:\n');
fprintf(fid, '  G4: %.1f Hz (+%.0f cents)\n', config.target_g4, cents_g4);
fprintf(fid, '  F3: %.1f Hz (%.0f cents)\n\n', config.target_f3, cents_f3);


fprintf('Pitch shifts: G4 = %.0f cents, F3 = %.0f cents\n', cents_g4, cents_f3);

% Parameter sweep values from config
sweep_config = project_config('exp1_sweep');
ola_synHop_values = sweep_config.ola_synHop_values;
ola_winLen_values = sweep_config.ola_winLen_values;
ola_winType_values = sweep_config.ola_winType_values;
wsola_tolerance_values = sweep_config.wsola_tolerance_values;

fprintf('\nStep 2a: OLA parameter sweep...\n');
fprintf(fid, 'OLA PARAMETER SWEEP:\n');
fprintf(fid, '==================\n\n');

best_ola_error = inf;
best_ola_params = [];
best_ola_signals = [];

% Sweep OLA parameters
for synHop = ola_synHop_values
    for winLen = ola_winLen_values
        for winType = ola_winType_values
            
            fprintf('Testing OLA: synHop=%d, winLen=%d, winType=%d\n', synHop, winLen, winType);
            fprintf(fid, 'Testing OLA: synHop=%d, winLen=%d, winType=%d\n', synHop, winLen, winType);
            
            % Set up OLA parameters
            param_ola.fsAudio = fs;
            param_ola.algTSM = @wsolaTSM;
            param_ola.tolerance = 0;
            param_ola.synHop = synHop;
            param_ola.win = win(winLen, winType);
            
            try
                % Test both pitch shifts
                y_ola_g4_test = pitchShiftViaTSM(x_original, cents_g4, param_ola);
                y_ola_f3_test = pitchShiftViaTSM(x_original, cents_f3, param_ola);
                
                % Analyze frequency accuracy
                [~, ~, peaks_g4] = analyze_fft(y_ola_g4_test, fs, false);
                [~, ~, peaks_f3] = analyze_fft(y_ola_f3_test, fs, false);
                
                % Save plots for this parameter combination
                param_str = sprintf('synHop=%d winLen=%d winType=%d', synHop, winLen, winType);
                
                % Save amplitude plots
                signals_g4_sweep = {x_original, y_ola_g4_test};
                titles_g4_sweep = {'Original (C4+G5)', sprintf('OLA G4 (%s)', param_str)};
                plot_amplitude_vs_time(signals_g4_sweep, t, titles_g4_sweep, [0, 0.02]);
                saveas(gcf, fullfile(plots_dir, sprintf('sweep_ola_g4_%s.png', param_str)));
                
                signals_f3_sweep = {x_original, y_ola_f3_test};
                titles_f3_sweep = {'Original (C4+G5)', sprintf('OLA F3 (%s)', param_str)};
                plot_amplitude_vs_time(signals_f3_sweep, t, titles_f3_sweep, [0, 0.02]);
                saveas(gcf, fullfile(plots_dir, sprintf('sweep_ola_f3_%s.png', param_str)));
                
                % Save FFT plots
                [~, ~, ~] = analyze_fft(y_ola_g4_test, fs, false);
                title(sprintf('OLA G4 FFT (%s)', param_str));
                saveas(gcf, fullfile(plots_dir, sprintf('sweep_ola_g4_fft_%s.png', param_str)));
                
                [~, ~, ~] = analyze_fft(y_ola_f3_test, fs, false);
                title(sprintf('OLA F3 FFT (%s)', param_str));
                saveas(gcf, fullfile(plots_dir, sprintf('sweep_ola_f3_fft_%s.png', param_str)));
                
                % Calculate combined error
                if ~isempty(peaks_g4.frequencies) && ~isempty(peaks_f3.frequencies)
                    error_g4 = abs(peaks_g4.frequencies(1) - config.target_g4);
                    error_f3 = abs(peaks_f3.frequencies(1) - config.target_f3);
                    total_error = error_g4 + error_f3;
                    
                    fprintf('  -> G4: %.1f Hz (target: %.1f Hz, error: %.1f Hz)\n', peaks_g4.frequencies(1), config.target_g4, error_g4);
                    fprintf('  -> F3: %.1f Hz (target: %.1f Hz, error: %.1f Hz)\n', peaks_f3.frequencies(1), config.target_f3, error_f3);
                    fprintf('  -> Total error: %.1f Hz\n', total_error);
                    
                    fprintf(fid, '  Results: G4=%.1f Hz (error: %.1f), F3=%.1f Hz (error: %.1f), Total=%.1f Hz\n', ...
                        peaks_g4.frequencies(1), error_g4, peaks_f3.frequencies(1), error_f3, total_error);
                    
                    % Check if this is the best so far
                    if total_error < best_ola_error
                        best_ola_error = total_error;
                        best_ola_params = struct('synHop', synHop, 'winLen', winLen, 'winType', winType);
                        best_ola_signals = struct('g4', y_ola_g4_test, 'f3', y_ola_f3_test); % Save the signals
                        fprintf('  -> NEW BEST!\n');
                        fprintf(fid, '  -> NEW BEST!\n');
                    end
                else
                    fprintf('  -> Error: No frequency peaks detected\n');
                    fprintf(fid, '  -> Error: No frequency peaks detected\n');
                end
            catch ME
                fprintf('  -> Error: %s\n', ME.message);
                fprintf(fid, '  -> Error: %s\n', ME.message);
            end
            fprintf('\n');
            fprintf(fid, '\n');
        end
    end
end

fprintf('Best OLA parameters: synHop=%d, winLen=%d, winType=%d (error=%.1f Hz)\n', ...
    best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType, best_ola_error);

fprintf(fid, '\nBEST OLA PARAMETERS:\n');
fprintf(fid, 'synHop=%d, winLen=%d, winType=%d\n', best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType);
fprintf(fid, 'Total error: %.1f Hz\n\n', best_ola_error);

fprintf('\nStep 2b: WSOLA tolerance sweep (using best OLA parameters)...\n');
fprintf(fid, 'WSOLA TOLERANCE SWEEP:\n');
fprintf(fid, '=====================\n');
fprintf(fid, 'Using best OLA parameters: synHop=%d, winLen=%d, winType=%d\n\n', ...
    best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType);

best_wsola_error = inf;
best_tolerance = [];
best_wsola_signals = [];

% Use best OLA parameters as base for WSOLA
for tolerance = wsola_tolerance_values
    
    fprintf('Testing WSOLA: tolerance=%d (using OLA params: synHop=%d, winLen=%d, winType=%d)\n', ...
        tolerance, best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType);
    fprintf(fid, 'Testing WSOLA: tolerance=%d\n', tolerance);
    
    % Set up WSOLA parameters (using best OLA params)
    param_wsola.fsAudio = fs;
    param_wsola.algTSM = @wsolaTSM;
    param_wsola.tolerance = tolerance;
    param_wsola.synHop = best_ola_params.synHop;
    param_wsola.win = win(best_ola_params.winLen, best_ola_params.winType);
    
    try
        % Test both pitch shifts
        y_wsola_g4_test = pitchShiftViaTSM(x_original, cents_g4, param_wsola);
        y_wsola_f3_test = pitchShiftViaTSM(x_original, cents_f3, param_wsola);
        
        % Analyze frequency accuracy
        [~, ~, peaks_g4] = analyze_fft(y_wsola_g4_test, fs, false);
        [~, ~, peaks_f3] = analyze_fft(y_wsola_f3_test, fs, false);
        
        % Save plots for this tolerance value
        tolerance_str = sprintf('tolerance=%d', tolerance);
        
        % Save amplitude plots
        signals_g4_sweep = {x_original, y_wsola_g4_test};
        titles_g4_sweep = {'Original (C4+G5)', sprintf('WSOLA G4 (%s)', tolerance_str)};
        plot_amplitude_vs_time(signals_g4_sweep, t, titles_g4_sweep, [0, 0.02]);
        saveas(gcf, fullfile(plots_dir, sprintf('sweep_wsola_g4_%s.png', tolerance_str)));
        
        signals_f3_sweep = {x_original, y_wsola_f3_test};
        titles_f3_sweep = {'Original (C4+G5)', sprintf('WSOLA F3 (%s)', tolerance_str)};
        plot_amplitude_vs_time(signals_f3_sweep, t, titles_f3_sweep, [0, 0.02]);
        saveas(gcf, fullfile(plots_dir, sprintf('sweep_wsola_f3_%s.png', tolerance_str)));
        
        % Save FFT plots
        [~, ~, ~] = analyze_fft(y_wsola_g4_test, fs, false);
        title(sprintf('WSOLA G4 FFT (%s)', tolerance_str));
        saveas(gcf, fullfile(plots_dir, sprintf('sweep_wsola_g4_fft_%s.png', tolerance_str)));
        
        [~, ~, ~] = analyze_fft(y_wsola_f3_test, fs, false);
        title(sprintf('WSOLA F3 FFT (%s)', tolerance_str));
        saveas(gcf, fullfile(plots_dir, sprintf('sweep_wsola_f3_fft_%s.png', tolerance_str)));
        
        % Calculate combined error
        if ~isempty(peaks_g4.frequencies) && ~isempty(peaks_f3.frequencies)
            error_g4 = abs(peaks_g4.frequencies(1) - config.target_g4);
            error_f3 = abs(peaks_f3.frequencies(1) - config.target_f3);
            total_error = error_g4 + error_f3;
            
            fprintf('  -> G4: %.1f Hz (target: %.1f Hz, error: %.1f Hz)\n', peaks_g4.frequencies(1), config.target_g4, error_g4);
            fprintf('  -> F3: %.1f Hz (target: %.1f Hz, error: %.1f Hz)\n', peaks_f3.frequencies(1), config.target_f3, error_f3);
            fprintf('  -> Total error: %.1f Hz\n', total_error);
            
            fprintf(fid, '  Results: G4=%.1f Hz (error: %.1f), F3=%.1f Hz (error: %.1f), Total=%.1f Hz\n', ...
                peaks_g4.frequencies(1), error_g4, peaks_f3.frequencies(1), error_f3, total_error);
            
            % Check if this is the best so far
            if total_error < best_wsola_error
                best_wsola_error = total_error;
                best_tolerance = tolerance;
                best_wsola_signals = struct('g4', y_wsola_g4_test, 'f3', y_wsola_f3_test); % Save the signals
                fprintf('  -> NEW BEST!\n');
                fprintf(fid, '  -> NEW BEST!\n');
            end
        else
            fprintf('  -> Error: No frequency peaks detected\n');
            fprintf(fid, '  -> Error: No frequency peaks detected\n');
        end
    catch ME
        fprintf('  -> Error: %s\n', ME.message);
        fprintf(fid, '  -> Error: %s\n', ME.message);
    end
    fprintf('\n');
    fprintf(fid, '\n');
end

fprintf('Best WSOLA tolerance: %d (error=%.1f Hz)\n', best_tolerance, best_wsola_error);

fprintf(fid, '\nBEST WSOLA PARAMETERS:\n');
fprintf(fid, 'tolerance=%d\n', best_tolerance);
fprintf(fid, 'Total error: %.1f Hz\n\n', best_wsola_error);

% Summary comparison
fprintf(fid, 'FINAL COMPARISON:\n');
fprintf(fid, '================\n');
fprintf(fid, 'OLA (optimal):    synHop=%d, winLen=%d, winType=%d, tolerance=0    -> Error: %.1f Hz\n', ...
    best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType, best_ola_error);
fprintf(fid, 'WSOLA (optimal):  synHop=%d, winLen=%d, winType=%d, tolerance=%d -> Error: %.1f Hz\n', ...
    best_ola_params.synHop, best_ola_params.winLen, best_ola_params.winType, best_tolerance, best_wsola_error);
fprintf(fid, '\nImprovement: %.1f Hz (%.1f%% better)\n', ...
    best_ola_error - best_wsola_error, ((best_ola_error - best_wsola_error) / best_ola_error) * 100);

% Close log file
fclose(fid);
fprintf('Parameter sweep log saved to: %s\n', log_file);

% Step 3: Use optimal parameters (signals already computed during sweep)
fprintf('\nStep 3: Using optimal parameters found during sweep...\n');

% Use the signals from the best parameter combinations (already computed)
y_ola_g4 = best_ola_signals.g4;
y_ola_f3 = best_ola_signals.f3;
y_wsola_g4 = best_wsola_signals.g4;
y_wsola_f3 = best_wsola_signals.f3;

% Step 4: Plot amplitude vs time comparisons
fprintf('\nStep 4: Plotting amplitude vs time comparisons...\n');

% Plot comparison for G4 pitch shift
signals_g4 = {x_original, y_ola_g4, y_wsola_g4};
titles_g4 = {'Original (C4+G5)', sprintf('OLA → G4 (%.0f cents)', cents_g4), sprintf('WSOLA → G4 (%.0f cents)', cents_g4)};
plot_amplitude_vs_time(signals_g4, t, titles_g4, [0, 0.02]);
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_G4.png'));

% Plot comparison for F3 pitch shift
signals_f3 = {x_original, y_ola_f3, y_wsola_f3};
titles_f3 = {'Original (C4+G5)', sprintf('OLA → F3 (%.0f cents)', cents_f3), sprintf('WSOLA → F3 (%.0f cents)', cents_f3)};
plot_amplitude_vs_time(signals_f3, t, titles_f3, [0, 0.02]);
saveas(gcf, fullfile(plots_dir, 'amplitude_vs_time_F3.png'));

% Step 5: FFT Analysis
fprintf('\nStep 5: Analyzing frequency content...\n');

% Analyze original signal
fprintf('\n--- ORIGINAL SIGNAL (C4 + G5) ---\n');
[~, ~, peaks_orig] = analyze_fft(x_original, fs, false);
title('Original Signal FFT (C4+G5)');
saveas(gcf, fullfile(plots_dir, 'fft_original_signal.png'));

% Analyze OLA results
fprintf('\n--- OLA G4 PITCH SHIFT (%.0f cents) ---\n', cents_g4);
[~, ~, peaks_ola_g4] = analyze_fft(y_ola_g4, fs, false);
title(sprintf('OLA G4 Pitch Shift FFT (%.0f cents)', cents_g4));
saveas(gcf, fullfile(plots_dir, 'fft_ola_g4.png'));

fprintf('\n--- OLA F3 PITCH SHIFT (%.0f cents) ---\n', cents_f3);
[~, ~, peaks_ola_f3] = analyze_fft(y_ola_f3, fs, false);
title(sprintf('OLA F3 Pitch Shift FFT (%.0f cents)', cents_f3));
saveas(gcf, fullfile(plots_dir, 'fft_ola_f3.png'));

% Analyze WSOLA results
fprintf('\n--- WSOLA G4 PITCH SHIFT (%.0f cents) ---\n', cents_g4);
[~, ~, peaks_wsola_g4] = analyze_fft(y_wsola_g4, fs, false);
title(sprintf('WSOLA G4 Pitch Shift FFT (%.0f cents)', cents_g4));
saveas(gcf, fullfile(plots_dir, 'fft_wsola_g4.png'));

fprintf('\n--- WSOLA F3 PITCH SHIFT (%.0f cents) ---\n', cents_f3);
[~, ~, peaks_wsola_f3] = analyze_fft(y_wsola_f3, fs, false);
title(sprintf('WSOLA F3 Pitch Shift FFT (%.0f cents)', cents_f3));
saveas(gcf, fullfile(plots_dir, 'fft_wsola_f3.png'));

fprintf('FFT plots saved to: %s\n', plots_dir);

% Step 6: Save audio files
fprintf('\nStep 6: Saving audio files...\n');

output_dir = fullfile('outputs', 'experiment1', 'audio');

% Save all signals
audiowrite(fullfile(output_dir, 'original_C4_G5.wav'), x_original, fs);
audiowrite(fullfile(output_dir, sprintf('OLA_G4_%.0f_cents.wav', cents_g4)), y_ola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('OLA_F3_%.0f_cents.wav', cents_f3)), y_ola_f3, fs);
audiowrite(fullfile(output_dir, sprintf('WSOLA_G4_%.0f_cents.wav', cents_g4)), y_wsola_g4, fs);
audiowrite(fullfile(output_dir, sprintf('WSOLA_F3_%.0f_cents.wav', cents_f3)), y_wsola_f3, fs);

fprintf('Audio files saved to: %s\n', output_dir);

% Step 7: Store results
results.optimal_params.ola = best_ola_params;
results.optimal_params.wsola_tolerance = best_tolerance;
results.optimal_params.error_ola = best_ola_error;
results.optimal_params.error_wsola = best_wsola_error;
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

