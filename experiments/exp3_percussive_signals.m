function results = exp3_percussive_signals()
% EXP3_PERCUSSIVE_SIGNALS - Experiment 3: Percussive Signal Processing (Bongo)
%
% This experiment analyzes how WSOLA and OLA handle percussive signals,
% specifically examining transient doubling/stuttering artifacts that are
% characteristic of time-domain TSM methods when processing signals with sharp onsets.
%
% Process:
%   1. Loads Bongo_ORIG.wav (fs = 22050 Hz)
%   2. Tests WSOLA with different tolerance values (0, 128, 256, 512, 1024, 2048)
%   3. Tests three stretching factors: α = 0.5 (compression), α = 1.0 (identity), α = 2.0 (stretching)
%   4. Compares WSOLA results against OLA baseline (tolerance = 0)
%   5. Demonstrates transient doubling/stuttering (α > 1) and transient skipping (α < 1)
%   6. Generates amplitude vs time plots focusing on transient behavior
%   7. Saves audio files showing different artifact types
%
% Focus:
%   - Transient doubling/stuttering artifacts with WSOLA
%   - Comparison showing OLA actually handles percussive content better
%   - Impact of tolerance parameter on artifact severity
%
% Output:
%   results - Structure containing:
%             .config - Experiment configuration
%             .original_signal - Input bongo signal
%             .processed_signals - All WSOLA/OLA results organized by alpha and tolerance
%             .analysis_summary - Key findings about transient preservation

fprintf('=== Experiment 3: Percussive Signal Processing (Bongo) ===\n\n');

% Step 1: Load configuration and input signal
fprintf('Step 1: Loading bongo signal...\n');
config = project_config('exp3');

% Try to load the bongo file
audio_dir = 'audio_inputs';
bongo_path = fullfile(audio_dir, config.input_file);

if ~exist(bongo_path, 'file')
    error('Bongo file not found: %s\nPlease ensure %s exists in the audio_inputs directory.', ...
          bongo_path, config.input_file);
end

[x_original, fs] = audioread(bongo_path);

% Verify sampling rate
if fs ~= config.expected_fs
    fprintf('Warning: Expected fs=%d Hz, got fs=%d Hz. Resampling...\n', config.expected_fs, fs);
    x_original = resample(x_original, config.expected_fs, fs);
    fs = config.expected_fs;
end

% Convert to mono if stereo
if size(x_original, 2) > 1
    x_original = mean(x_original, 2);
    fprintf('Converted stereo to mono\n');
end

% Limit duration for analysis
max_samples = round(config.analysis_duration * fs);
if length(x_original) > max_samples
    x_original = x_original(1:max_samples);
    fprintf('Truncated to %.1f seconds for analysis\n', config.analysis_duration);
end

fprintf('Loaded: %s (%.2f seconds, %d Hz)\n', config.input_file, length(x_original)/fs, fs);

% Create output directories
plots_dir = fullfile('outputs', 'experiment3', 'plots');
audio_dir_out = fullfile('outputs', 'experiment3', 'audio');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
if ~exist(audio_dir_out, 'dir'), mkdir(audio_dir_out); end

% Initialize results structure
results = struct();
results.config = config;
results.fs = fs;
results.original_signal = x_original;
results.processed_signals = struct();

% Step 2: Test all combinations of alpha and tolerance
fprintf('\nStep 2: Testing WSOLA with different parameters...\n');

% Create log file
log_file = fullfile('outputs', 'experiment3', 'processing_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 3: Percussive Signal Processing Log ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));

for alpha_idx = 1:length(config.alpha_values)
    alpha = config.alpha_values(alpha_idx);
    
    fprintf('\n--- Testing Alpha = %.1f ---\n', alpha);
    fprintf(fid, '\nAlpha = %.1f:\n', alpha);
    fprintf(fid, '================\n');
    
    % Initialize storage for this alpha
    alpha_field = sprintf('alpha_%.1f', alpha);
    alpha_field = strrep(alpha_field, '.', '_');
    results.processed_signals.(alpha_field) = struct();
    
    for tol_idx = 1:length(config.tolerance_values)
        tolerance = config.tolerance_values(tol_idx);
        
        fprintf('  Testing tolerance = %d samples... ', tolerance);
        fprintf(fid, '  Tolerance = %d samples: ', tolerance);
        
        % Set up WSOLA parameters
        param_wsola = struct();
        param_wsola.tolerance = tolerance;
        param_wsola.synHop = config.syn_hop;
        param_wsola.win = win(config.frame_size, config.window_beta);
        
        try
            % Apply WSOLA
            y_wsola = wsolaTSM(x_original, alpha, param_wsola);
            
            % Store result
            tol_field = sprintf('tol_%d', tolerance);
            results.processed_signals.(alpha_field).(tol_field) = y_wsola;
            
            % Log success
            fprintf('SUCCESS (length: %d samples)\n', length(y_wsola));
            fprintf(fid, 'SUCCESS (output length: %d samples)\n', length(y_wsola));
            
            % Save audio file immediately
            save_audio_file(y_wsola, fs, alpha, tolerance, audio_dir_out);
            
            % Generate individual plot immediately
            create_individual_plot(x_original, y_wsola, fs, alpha, tolerance, config.plot_time_range, plots_dir);
            
        catch ME
            fprintf('FAILED: %s\n', ME.message);
            fprintf(fid, 'FAILED: %s\n', ME.message);
            results.processed_signals.(alpha_field).(tol_field) = [];
        end
    end
end

% Close log file
fclose(fid);
fprintf('Processing log saved to: %s\n', log_file);

% Step 3: Generate comparison plots
fprintf('\nStep 3: Generating comparison plots...\n');

% Create tolerance comparison plots for each alpha
for alpha_idx = 1:length(config.alpha_values)
    alpha = config.alpha_values(alpha_idx);
    alpha_field = sprintf('alpha_%.1f', alpha);
    alpha_field = strrep(alpha_field, '.', '_');
    
    if isfield(results.processed_signals, alpha_field)
        create_tolerance_comparison(x_original, results.processed_signals.(alpha_field), ...
                                  fs, alpha, config, plots_dir);
    end
end

% Create overall summary plot
fprintf('\nStep 4: Creating summary analysis plot...\n');
create_ola_wsola_summary(x_original, results.processed_signals, fs, config, plots_dir);

% Step 5: Save original signal for reference
fprintf('\nStep 5: Saving reference files...\n');
audiowrite(fullfile(audio_dir_out, 'original_bongo.wav'), x_original, fs);

% Generate analysis summary
results.analysis_summary = create_analysis_summary(results.processed_signals, config);

fprintf('\nExperiment 3 completed successfully!\n');
fprintf('Audio files saved to: %s\n', audio_dir_out);
fprintf('Plots saved to: %s\n', plots_dir);
fprintf('Results stored in output structure.\n');

end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function save_audio_file(y_processed, fs, alpha, tolerance, audio_dir)
% Save processed audio file with descriptive name and normalization

% Normalize to prevent clipping
max_val = max(abs(y_processed));
if max_val > 0.95
    y_processed = y_processed * (0.95 / max_val);
end

if tolerance == 0
    filename = sprintf('bongo_OLA_alpha_%.1f_tol_%d.wav', alpha, tolerance);
else
    filename = sprintf('bongo_WSOLA_alpha_%.1f_tol_%d.wav', alpha, tolerance);
end

% Replace decimal points in numbers but preserve the .wav extension
filename = strrep(filename, '.wav', '_TEMP_WAV');
filename = strrep(filename, '.', '_');
filename = strrep(filename, '_TEMP_WAV', '.wav');

audiowrite(fullfile(audio_dir, filename), y_processed, fs);

end

function create_individual_plot(x_original, y_processed, fs, alpha, tolerance, time_range, plots_dir)
% Generate and save individual comparison plot

% Create time vectors
t_orig = (0:length(x_original)-1) / fs;
t_proc = (0:length(y_processed)-1) / fs;

% Create figure
figure('Position', [100, 100, 800, 400], 'Visible', 'off');

% Plot original
subplot(2, 1, 1);
plot(t_orig, x_original, 'b-', 'LineWidth', 1);
grid on;
title('Original Bongo Signal', 'FontSize', 12);
xlabel('Time (s)');
ylabel('Amplitude');
xlim(time_range);
ylim([-1.1, 1.1]);

% Plot processed
subplot(2, 1, 2);
plot(t_proc, y_processed, 'r-', 'LineWidth', 1);
grid on;
if tolerance == 0
    title(sprintf('OLA (α=%.1f, tolerance=%d)', alpha, tolerance), 'FontSize', 12);
else
    title(sprintf('WSOLA (α=%.1f, tolerance=%d)', alpha, tolerance), 'FontSize', 12);
end
xlabel('Time (s)');
ylabel('Amplitude');
xlim(time_range * alpha);
ylim([-1.1, 1.1]);

% Add overall title
sgtitle(sprintf('Percussive Signal Processing: α=%.1f, tolerance=%d', alpha, tolerance), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Save plot
plot_filename = sprintf('individual_alpha_%.1f_tol_%d.png', alpha, tolerance);
plot_filename = strrep(plot_filename, '.png', '_TEMP_PNG');
plot_filename = strrep(plot_filename, '.', '_');
plot_filename = strrep(plot_filename, '_TEMP_PNG', '.png');
saveas(gcf, fullfile(plots_dir, plot_filename));
close(gcf);

end

function create_tolerance_comparison(x_original, alpha_results, fs, alpha, config, plots_dir)
% Create comparison plot for all tolerance values at specific alpha

% Collect successful results
tolerance_values = config.tolerance_values;
signals = {x_original};
titles = {'Original'};

for i = 1:length(tolerance_values)
    tol_field = sprintf('tol_%d', tolerance_values(i));
    if isfield(alpha_results, tol_field) && ~isempty(alpha_results.(tol_field))
        signals{end+1} = alpha_results.(tol_field);
        if tolerance_values(i) == 0
            titles{end+1} = sprintf('OLA (tol=%d)', tolerance_values(i));
        else
            titles{end+1} = sprintf('WSOLA (tol=%d)', tolerance_values(i));
        end
    end
end

if length(signals) > 1
    % Find the longest signal for time vector
    max_length = max(cellfun(@length, signals));
    t_common = (0:max_length-1) / fs;
    
    % Pad all signals to same length
    padded_signals = cell(size(signals));
    for i = 1:length(signals)
        signal = signals{i};
        padded_signals{i} = [signal; zeros(max_length - length(signal), 1)];
    end
    
    % Create plot using the utility function
    fig = plot_amplitude_vs_time(padded_signals, t_common, titles, config.plot_time_range, false);
    sgtitle(sprintf('Tolerance Comparison: α=%.1f', alpha), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save plot
    filename = sprintf('alpha_%.1f_tolerance_comparison.png', alpha);
    filename = strrep(filename, '.png', '_TEMP_PNG');
    filename = strrep(filename, '.', '_');
    filename = strrep(filename, '_TEMP_PNG', '.png');
    saveas(gcf, fullfile(plots_dir, filename));
    close(gcf);
end

end

function create_ola_wsola_summary(x_original, processed_signals, fs, config, plots_dir)
% Create overall summary comparing OLA vs best WSOLA for each alpha

figure('Position', [100, 100, 1200, 800], 'Visible', 'off');

alpha_values = config.alpha_values;
plot_count = 0;

for alpha_idx = 1:length(alpha_values)
    alpha = alpha_values(alpha_idx);
    alpha_field = sprintf('alpha_%.1f', alpha);
    alpha_field = strrep(alpha_field, '.', '_');
    
    if isfield(processed_signals, alpha_field)
        alpha_results = processed_signals.(alpha_field);
        
        % Get OLA result (tolerance = 0)
        if isfield(alpha_results, 'tol_0') && ~isempty(alpha_results.tol_0)
            y_ola = alpha_results.tol_0;
            
            % Get a representative WSOLA result (try tolerance = 512 first)
            y_wsola = [];
            wsola_tol = 0;
            for test_tol = [512, 256, 1024, 128, 2048]
                tol_field = sprintf('tol_%d', test_tol);
                if isfield(alpha_results, tol_field) && ~isempty(alpha_results.(tol_field))
                    y_wsola = alpha_results.(tol_field);
                    wsola_tol = test_tol;
                    break;
                end
            end
            
            if ~isempty(y_wsola)
                plot_count = plot_count + 1;
                
                % Create time vectors
                t_orig = (0:length(x_original)-1) / fs;
                t_ola = (0:length(y_ola)-1) / fs;
                t_wsola = (0:length(y_wsola)-1) / fs;
                
                % Original signal
                subplot(length(alpha_values), 3, (plot_count-1)*3 + 1);
                plot(t_orig, x_original, 'k-', 'LineWidth', 1);
                grid on;
                title(sprintf('Original (α=%.1f)', alpha));
                xlabel('Time (s)');
                ylabel('Amplitude');
                xlim(config.plot_time_range);
                ylim([-1.1, 1.1]);
                
                % OLA result
                subplot(length(alpha_values), 3, (plot_count-1)*3 + 2);
                plot(t_ola, y_ola, 'b-', 'LineWidth', 1);
                grid on;
                title(sprintf('OLA (α=%.1f)', alpha));
                xlabel('Time (s)');
                ylabel('Amplitude');
                xlim(config.plot_time_range * alpha);
                ylim([-1.1, 1.1]);
                
                % WSOLA result
                subplot(length(alpha_values), 3, (plot_count-1)*3 + 3);
                plot(t_wsola, y_wsola, 'r-', 'LineWidth', 1);
                grid on;
                title(sprintf('WSOLA (α=%.1f, tol=%d)', alpha, wsola_tol));
                xlabel('Time (s)');
                ylabel('Amplitude');
                xlim(config.plot_time_range * alpha);
                ylim([-1.1, 1.1]);
            end
        end
    end
end

% Add overall title
sgtitle('Percussive Signal Processing: OLA vs WSOLA Summary', 'FontSize', 16, 'FontWeight', 'bold');

% Save plot
saveas(gcf, fullfile(plots_dir, 'summary_ola_vs_wsola.png'));
close(gcf);

end

function summary = create_analysis_summary(processed_signals, config)
% Generate analysis summary of the experiment results

summary = struct();
summary.alpha_values_tested = config.alpha_values;
summary.tolerance_values_tested = config.tolerance_values;
summary.successful_results = struct();

% Count successful results for each alpha
alpha_values = config.alpha_values;
for alpha_idx = 1:length(alpha_values)
    alpha = alpha_values(alpha_idx);
    alpha_field = sprintf('alpha_%.1f', alpha);
    alpha_field = strrep(alpha_field, '.', '_');
    
    if isfield(processed_signals, alpha_field)
        alpha_results = processed_signals.(alpha_field);
        
        successful_count = 0;
        tolerance_list = [];
        
        for tol_idx = 1:length(config.tolerance_values)
            tolerance = config.tolerance_values(tol_idx);
            tol_field = sprintf('tol_%d', tolerance);
            
            if isfield(alpha_results, tol_field) && ~isempty(alpha_results.(tol_field))
                successful_count = successful_count + 1;
                tolerance_list(end+1) = tolerance;
            end
        end
        
        summary.successful_results.(alpha_field) = struct();
        summary.successful_results.(alpha_field).count = successful_count;
        summary.successful_results.(alpha_field).tolerance_values = tolerance_list;
    end
end

end