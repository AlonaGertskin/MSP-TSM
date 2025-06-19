function results = exp5_polyphonic_music()
% EXP5_POLYPHONIC_MUSIC - Experiment 5: Polyphonic Music Layers
%
% This experiment analyzes how OLA and WSOLA handle polyphonic music with
% different layers and combinations. Tests the effect of window sizes on
% different musical content types (harmonic vs percussive).
%
% Process:
%   1. Loads individual layers (bass, drums, pads, synth)
%   2. Loads layer combinations and full mix
%   3. Tests at 150% and 75% speeds with multiple window sizes
%   4. Compares OLA vs WSOLA performance across different content types
%   5. Analyzes layer-specific sensitivity to parameters
%
% Output:
%   results - Structure containing all processed signals and analysis

fprintf('=== Experiment 5: Polyphonic Music Layers ===\n\n');

% Load configuration
config = project_config('exp5');

% Create output directories
plots_dir = fullfile('outputs', 'experiment5', 'plots');
audio_dir = fullfile('outputs', 'experiment5', 'audio');

% Initialize results structure
results = struct();
results.config = config;
results.fs = config.expected_fs;
results.processed_signals = struct();

% Step 1: Load all audio files
fprintf('Step 1: Loading audio files...\n');

% Initialize audio storage
audio_files = struct();
audio_files.individual = struct();
audio_files.combinations = struct();
audio_files.full_mix = struct();

% Load individual layers
fprintf('  Loading individual layers...\n');
for i = 1:length(config.individual_layers)
    filename = config.individual_layers{i};
    filepath = fullfile('audio_inputs', filename);
    
    if ~exist(filepath, 'file')
        warning('File not found: %s', filepath);
        continue;
    end
    
    [audio_data, fs] = audioread(filepath);
    
    % Handle sampling rate consistency
    if fs ~= config.expected_fs
        fprintf('    Resampling %s from %d Hz to %d Hz\n', filename, fs, config.expected_fs);
        audio_data = resample(audio_data, config.expected_fs, fs);
        fs = config.expected_fs;
    end
    
    % Convert to mono if needed
    if size(audio_data, 2) > 1
        audio_data = mean(audio_data, 2);
        fprintf('    Converted %s to mono\n', filename);
    end
    
    % Skip initial silence
    if length(audio_data) > round(config.start_offset_seconds * fs)
        start_samples = round(config.start_offset_seconds * fs);
        audio_data = audio_data(start_samples+1:end);
    end

    % Limit duration
    max_samples = round(config.analysis_duration * fs);
    if length(audio_data) > max_samples
        audio_data = audio_data(1:max_samples);
    end
    
    % Store with clean field name
    field_name = strrep(strrep(filename, '.wav', ''), '_only', '');
    audio_files.individual.(field_name) = audio_data;
    fprintf('    Loaded: %s (%.1f seconds)\n', filename, length(audio_data)/fs);
end

% Load layer combinations
fprintf('  Loading layer combinations...\n');
for i = 1:length(config.layer_combinations)
    filename = config.layer_combinations{i};
    filepath = fullfile('audio_inputs', filename);
    
    if ~exist(filepath, 'file')
        warning('File not found: %s', filepath);
        continue;
    end
    
    [audio_data, fs] = audioread(filepath);
    
    % Same processing as individual layers (sampling rate, mono, duration)
    if fs ~= config.expected_fs
        audio_data = resample(audio_data, config.expected_fs, fs);
        fs = config.expected_fs;
    end
    
    if size(audio_data, 2) > 1
        audio_data = mean(audio_data, 2);
    end
    
    max_samples = round(config.analysis_duration * fs);
    if length(audio_data) > max_samples
        audio_data = audio_data(1:max_samples);
    end
    
    % Store with clean field name
    field_name = strrep(filename, '.wav', '');
    audio_files.combinations.(field_name) = audio_data;
    fprintf('    Loaded: %s (%.1f seconds)\n', filename, length(audio_data)/fs);
end

% Load full mix
fprintf('  Loading full mix...\n');
full_mix_path = fullfile('audio_inputs', config.full_mix_file);
if exist(full_mix_path, 'file')
    [audio_data, fs] = audioread(full_mix_path);
    
    % Same processing
    if fs ~= config.expected_fs
        audio_data = resample(audio_data, config.expected_fs, fs);
    end
    
    if size(audio_data, 2) > 1
        audio_data = mean(audio_data, 2);
    end
    
    max_samples = round(config.analysis_duration * fs);
    if length(audio_data) > max_samples
        audio_data = audio_data(1:max_samples);
    end
    
    field_name = strrep(config.full_mix_file, '.wav', '');
    audio_files.full_mix.(field_name) = audio_data;
    fprintf('    Loaded: %s (%.1f seconds)\n', config.full_mix_file, length(audio_data)/fs);
else
    warning('Full mix file not found: %s', full_mix_path);
end

fprintf('Audio loading completed.\n\n');

% Step 2: Setup processing parameters
fprintf('Step 2: Setting up processing parameters...\n');

% Convert window sizes from milliseconds to samples
window_sizes_samples = round(config.window_sizes_ms * config.expected_fs / 1000);
fprintf('  Window sizes: ');
for i = 1:length(config.window_sizes_ms)
    fprintf('%dms (%d samples) ', config.window_sizes_ms(i), window_sizes_samples(i));
end
fprintf('\n');

% Setup algorithm configurations
algorithms = struct();

% OLA configuration (tolerance = 0)
algorithms.OLA.name = 'OLA';
algorithms.OLA.tolerance = 0;
algorithms.OLA.synHop = config.syn_hop_size;

% WSOLA configuration  
algorithms.WSOLA.name = 'WSOLA';
algorithms.WSOLA.tolerance = config.tolerance;
algorithms.WSOLA.synHop = config.syn_hop_size;

% Window function (consistent for both algorithms)
win_function = win(config.frame_size, config.window_beta);  % Hann window

% Display processing matrix
fprintf('  Processing matrix:\n');
fprintf('    Algorithms: %s, %s\n', algorithms.OLA.name, algorithms.WSOLA.name);
fprintf('    Speeds: ');
for alpha = config.alpha_values
    if alpha > 1
        fprintf('%.0f%% speed (α=%.2f) ', alpha*100, alpha);
    else
        fprintf('%.0f%% speed (α=%.2f) ', alpha*100, alpha);
    end
end
fprintf('\n');

% Count total processing tasks
num_individual = length(fieldnames(audio_files.individual));
num_combinations = length(fieldnames(audio_files.combinations));
num_full_mix = ~isempty(audio_files.full_mix);
total_audio_files = num_individual + num_combinations + num_full_mix;

total_tasks = total_audio_files * length(fieldnames(algorithms)) * ...
              length(config.alpha_values) * length(window_sizes_samples);

fprintf('    Total audio files: %d (individual: %d, combinations: %d, full mix: %d)\n', ...
        total_audio_files, num_individual, num_combinations, num_full_mix);
fprintf('    Total processing tasks: %d\n', total_tasks);

% Initialize processing counter
processing_counter = 0;

fprintf('Parameter setup completed.\n\n');

% Step 3: Main processing loop
fprintf('Step 3: Processing audio files...\n');

% Create log file for detailed results
log_file = fullfile('outputs', 'experiment5', 'processing_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 5: Polyphonic Music Processing Log ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));

% Define processing order and categories
audio_categories = {'individual', 'combinations', 'full_mix'};
algorithm_names = fieldnames(algorithms);



% Main processing loops
for cat_idx = 1:length(audio_categories)
    category = audio_categories{cat_idx};
    fprintf('\n--- Processing %s files ---\n', category);
    fprintf(fid, '\n=== %s FILES ===\n', upper(category));
    
    % Get audio files for this category (same for all categories now)
    file_list = fieldnames(audio_files.(category));
    if isempty(file_list)
        fprintf('  No files in %s category\n', category);
        continue;
    end
    audio_data_list = struct2cell(audio_files.(category));
    
    % Process each file in this category
    for file_idx = 1:length(file_list)
        filename = file_list{file_idx};
        audio_data = audio_data_list{file_idx};
        
        fprintf('  Processing: %s\n', filename);
        fprintf(fid, '\nFile: %s\n', filename);
        fprintf(fid, 'Duration: %.2f seconds (%d samples)\n', length(audio_data)/config.expected_fs, length(audio_data));
        
        % Initialize storage for this file
        if ~isfield(results.processed_signals, category)
            results.processed_signals.(category) = struct();
        end
        results.processed_signals.(category).(filename) = struct();

        % Initialize storage for comparison plots
        comparison_data = struct();
        comparison_data.original = audio_data;
        comparison_data.results = struct();
        
        % Create plots directory for this file
        file_plots_dir = fullfile(plots_dir, filename);
        if ~exist(file_plots_dir, 'dir')
            mkdir(file_plots_dir);
        end
                % Algorithm loop
        for alg_idx = 1:length(algorithm_names)
            alg_name = algorithm_names{alg_idx};
            alg_config = algorithms.(alg_name);
            
            fprintf('    Algorithm: %s\n', alg_name);
            fprintf(fid, '  Algorithm: %s (tolerance=%d)\n', alg_name, alg_config.tolerance);
            
            % Initialize storage for this algorithm
            results.processed_signals.(category).(filename).(alg_name) = struct();
            
            % Speed loop
            for alpha_idx = 1:length(config.alpha_values)
                alpha = config.alpha_values(alpha_idx);
                speed_percent = alpha * 100;
                
                fprintf('      Speed: %.0f%% (α=%.2f)\n', speed_percent, alpha);
                fprintf(fid, '    Speed: %.0f%% (α=%.2f)\n', speed_percent, alpha);
                
                % Initialize storage for this speed
                speed_field = sprintf('speed_%.0f', speed_percent);
                results.processed_signals.(category).(filename).(alg_name).(speed_field) = struct();
                
                % Window size loop
                for win_idx = 1:length(window_sizes_samples)
                    window_samples = window_sizes_samples(win_idx);
                    window_ms = config.window_sizes_ms(win_idx);
                    
                    processing_counter = processing_counter + 1;
                    fprintf('        Window: %dms (%d/%d) ... ', window_ms, processing_counter, total_tasks);
                    fprintf(fid, '      Window: %dms (%d samples): ', window_ms, window_samples);
                    
                    % Create window field name
                    window_field = sprintf('win_%dms', window_ms);
                    
                    try
                        % Setup WSOLA parameters for this specific test
                        param_wsola = struct();
                        param_wsola.tolerance = alg_config.tolerance;
                        param_wsola.synHop = alg_config.synHop;
                        param_wsola.win = win(window_samples, config.window_beta);
                        
                        % Apply time-scale modification
                        processed_audio = wsolaTSM(audio_data, alpha, param_wsola);
                        
                        % Store the processed result
                        results.processed_signals.(category).(filename).(alg_name).(speed_field).(window_field) = processed_audio;
                        
                        % Create subdirectory for this audio file
                        file_audio_dir = fullfile(audio_dir, filename);
                        if ~exist(file_audio_dir, 'dir')
                            mkdir(file_audio_dir);
                        end
                        
                        % Generate descriptive filename for audio output
                        audio_filename = sprintf('%s_%s_speed%.0f_%s.wav', ...
                            alg_name, category, speed_percent, window_field);
                        
                        % Save to file-specific directory
                        output_path = fullfile(file_audio_dir, audio_filename);
                        
                        % Normalize to prevent clipping
                        if max(abs(processed_audio)) > 0
                            normalized_audio = processed_audio / max(abs(processed_audio)) * 0.95;
                        else
                            normalized_audio = processed_audio;
                        end
                        
                        audiowrite(output_path, normalized_audio, config.expected_fs);
                        
                        % Log success
                        fprintf('SUCCESS\n');
                        fprintf(fid, 'SUCCESS - Output: %d samples (%.2f sec)\n', ...
                            length(processed_audio), length(processed_audio)/config.expected_fs);
                        
                        % Generate amplitude plot for this specific case
                        create_individual_comparison_plot(audio_data, processed_audio, config.expected_fs, ...
                            filename, alg_name, alpha, window_ms, config.plot_time_range, plots_dir);

                        % Store result for comparison plots
                        if ~isfield(comparison_data.results, alg_name)
                            comparison_data.results.(alg_name) = struct();
                        end
                        if ~isfield(comparison_data.results.(alg_name), speed_field)
                            comparison_data.results.(alg_name).(speed_field) = struct();
                        end
                        comparison_data.results.(alg_name).(speed_field).(window_field) = processed_audio;
    
catch ME
    % Handle processing failures
    fprintf('FAILED: %s\n', ME.message);
    fprintf(fid, 'FAILED: %s\n', ME.message);
    
    % Store empty result to indicate failure
    results.processed_signals.(category).(filename).(alg_name).(speed_field).(window_field) = [];
end
                    
                end % window loop
            end % speed loop
        end % algorithm loop
% Generate comparison plots for this file
create_comparison_plots(comparison_data, filename, category, config, file_plots_dir);

        fprintf('    Completed: %s\n', filename);
        fprintf(fid, '  COMPLETED: %s\n', filename);
        
    end % file loop
end % category loop

% Close log file
fclose(fid);
fprintf('\nAll processing completed! Log saved to: %s\n', log_file);

end 

function create_individual_comparison_plot(original, processed, fs, filename, alg_name, alpha, window_ms, time_range, plots_dir)
% Create and save comparison plot for individual processing result

% Create time vectors
t_orig = (0:length(original)-1) / fs;
t_proc = (0:length(processed)-1) / fs;

% Create figure (hidden)
fig = figure('Position', [100, 100, 800, 400], 'Visible', 'off');

% Plot original
subplot(2, 1, 1);
plot(t_orig, original, 'b-', 'LineWidth', 1);
grid on;
title(sprintf('Original: %s', strrep(filename, '_', ' ')), 'FontSize', 10);
xlabel('Time (s)');
ylabel('Amplitude');
xlim(time_range);
ylim([-1.1, 1.1]);

% Plot processed
subplot(2, 1, 2);
plot(t_proc, processed, 'r-', 'LineWidth', 1);
grid on;
title(sprintf('%s: α=%.2f, Window=%dms', alg_name, alpha, window_ms), 'FontSize', 10);
xlabel('Time (s)');
ylabel('Amplitude');
xlim(time_range * alpha);  % Adjust time range for speed change
ylim([-1.1, 1.1]);

% Add overall title
sgtitle(sprintf('%s Processing: %s', alg_name, strrep(filename, '_', ' ')), 'FontSize', 12, 'FontWeight', 'bold');

% Create subdirectory for this audio file's plots
file_plots_dir = fullfile(plots_dir, filename);
if ~exist(file_plots_dir, 'dir')
    mkdir(file_plots_dir);
end

% Save plot to file-specific directory
plot_filename = sprintf('%s_speed%.0f_win%dms.png', alg_name, alpha*100, window_ms);
saveas(fig, fullfile(file_plots_dir, plot_filename));

end

function create_comparison_plots(comparison_data, filename, category, config, plots_dir)
% Create comparison plots showing different window sizes and algorithms

try
    original = comparison_data.original;
    results = comparison_data.results;
    fs = config.expected_fs;
    
    algorithm_names = fieldnames(results);
    
    % Create comparison for each speed
    for alpha = config.alpha_values
        speed_field = sprintf('speed_%.0f', alpha * 100);
        
        % Check if this speed has results
        has_results = false;
        for alg_idx = 1:length(algorithm_names)
            if isfield(results.(algorithm_names{alg_idx}), speed_field)
                has_results = true;
                break;
            end
        end
        
        if ~has_results, continue; end
        
        % Create window size comparison plot
        create_window_comparison(original, results, speed_field, alpha, filename, ...
                               config, fs, plots_dir);
        
        % Create algorithm comparison plot  
        create_algorithm_comparison(original, results, speed_field, alpha, filename, ...
                                  config, fs, plots_dir);
    end
    
catch ME
    warning('Comparison plot failed for %s: %s', filename, ME.message);
end

end

function create_window_comparison(original, results, speed_field, alpha, filename, config, fs, plots_dir)
% Compare different window sizes for each algorithm

algorithm_names = fieldnames(results);

for alg_idx = 1:length(algorithm_names)
    alg_name = algorithm_names{alg_idx};
    
    if ~isfield(results.(alg_name), speed_field), continue; end
    
    window_results = results.(alg_name).(speed_field);
    window_fields = fieldnames(window_results);
    
    if length(window_fields) < 2, continue; end % Need at least 2 windows to compare
    
    % Create figure
    fig = figure('Position', [100, 100, 1000, 600], 'Visible', 'off');
    
    num_plots = length(window_fields) + 1; % +1 for original
    
    % Plot original
    subplot(num_plots, 1, 1);
    t_orig = (0:length(original)-1) / fs;
    plot(t_orig, original, 'k-', 'LineWidth', 1);
    grid on;
    title(sprintf('Original: %s', strrep(filename, '_', ' ')), 'FontSize', 10);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim(config.plot_time_range);
    ylim([-1.1, 1.1]);
    
    % Plot each window size
    for w_idx = 1:length(window_fields)
        subplot(num_plots, 1, w_idx + 1);
        
        window_field = window_fields{w_idx};
        processed = window_results.(window_field);
        
        if isempty(processed), continue; end
        
        t_proc = (0:length(processed)-1) / fs;
        plot(t_proc, processed, 'b-', 'LineWidth', 1);
        grid on;
        
        % Extract window size from field name
        window_ms = str2double(strrep(strrep(window_field, 'win_', ''), 'ms', ''));
        title(sprintf('%s: %dms window (α=%.2f)', alg_name, window_ms, alpha), 'FontSize', 10);
        xlabel('Time (s)');
        ylabel('Amplitude');
        xlim(config.plot_time_range * alpha);
        ylim([-1.1, 1.1]);
    end
    
    sgtitle(sprintf('Window Size Comparison: %s (%s, %.0f%% speed)', ...
            strrep(filename, '_', ' '), alg_name, alpha*100), 'FontSize', 12, 'FontWeight', 'bold');
    
    % Save plot
    plot_filename = sprintf('%s_window_comparison_speed%.0f.png', alg_name, alpha*100);
    saveas(fig, fullfile(plots_dir, plot_filename));
    close(fig);
end

end

function create_algorithm_comparison(original, results, speed_field, alpha, filename, config, fs, plots_dir)
% Compare OLA vs WSOLA for each window size

window_sizes = {};
algorithm_names = fieldnames(results);

% Find common window sizes across algorithms
for alg_idx = 1:length(algorithm_names)
    alg_name = algorithm_names{alg_idx};
    if isfield(results.(alg_name), speed_field)
        if isempty(window_sizes)
            window_sizes = fieldnames(results.(alg_name).(speed_field));
        else
            % Keep only common window sizes
            current_windows = fieldnames(results.(alg_name).(speed_field));
            window_sizes = intersect(window_sizes, current_windows);
        end
    end
end

for w_idx = 1:length(window_sizes)
    window_field = window_sizes{w_idx};
    window_ms = str2double(strrep(strrep(window_field, 'win_', ''), 'ms', ''));
    
    % Create figure
    fig = figure('Position', [100, 100, 800, 600], 'Visible', 'off');
    
    num_plots = length(algorithm_names) + 1; % +1 for original
    
    % Plot original
    subplot(num_plots, 1, 1);
    t_orig = (0:length(original)-1) / fs;
    plot(t_orig, original, 'k-', 'LineWidth', 1);
    grid on;
    title(sprintf('Original: %s', strrep(filename, '_', ' ')), 'FontSize', 10);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim(config.plot_time_range);
    ylim([-1.1, 1.1]);
    
    % Plot each algorithm
    for alg_idx = 1:length(algorithm_names)
        subplot(num_plots, 1, alg_idx + 1);
        
        alg_name = algorithm_names{alg_idx};
        
        if isfield(results.(alg_name), speed_field) && ...
           isfield(results.(alg_name).(speed_field), window_field)
            
            processed = results.(alg_name).(speed_field).(window_field);
            
            if ~isempty(processed)
                t_proc = (0:length(processed)-1) / fs;
                plot(t_proc, processed, 'r-', 'LineWidth', 1);
                grid on;
                title(sprintf('%s (α=%.2f)', alg_name, alpha), 'FontSize', 10);
                xlabel('Time (s)');
                ylabel('Amplitude');
                xlim(config.plot_time_range * alpha);
                ylim([-1.1, 1.1]);
            else
                title(sprintf('%s: FAILED', alg_name), 'FontSize', 10);
            end
        else
            title(sprintf('%s: NO DATA', alg_name), 'FontSize', 10);
        end
    end
    
    sgtitle(sprintf('Algorithm Comparison: %s (%dms window, %.0f%% speed)', ...
            strrep(filename, '_', ' '), window_ms, alpha*100), 'FontSize', 12, 'FontWeight', 'bold');
    
    % Save plot
    plot_filename = sprintf('algorithm_comparison_%dms_speed%.0f.png', window_ms, alpha*100);
    saveas(fig, fullfile(plots_dir, plot_filename));
    close(fig);
end

end