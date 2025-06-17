function results = exp2_tempo_modification()
% EXP2_TEMPO_MODIFICATION - Experiment 2: Tempo Modification Analysis (Refactored)
%
% This refactored version reduces code duplication by extracting common
% functionality into helper functions.

fprintf('=== Experiment 2: Tempo Modification Analysis ===\n\n');

% Fixed sampling rate
fs = 44100;

% Step 1: Load configuration and generate test signal
fprintf('Step 1: Loading experiment configuration...\n');
config = project_config('exp2');
fprintf('Test signal generated: %.1f seconds, %.0f Hz + %.0f Hz\n', ...
    config.duration, config.freq1, config.freq2);

[x_original, t] = generate_periodic_signal(config.duration, config.freq1, config.freq2);

% Initialize results structure
results = struct();
results.config = config;
results.fs = fs;
results.original_signal = x_original;
results.time = t;
results.signals = struct();

% Step 2: Choose which parts to run
run_part_a = true;   % Hop Size Analysis (α = Hs/Ha)
run_part_b = true;  % Frame Size Effects
run_part_c = true;  % Window Type Comparison

if run_part_a
    fprintf('\nRunning Part A - Hop Size Analysis (Alpha Limits)...\n');
    results = run_experiment_part(x_original, t, fs, config, results, 'part_a');
end

if run_part_b
    fprintf('\nRunning Part B - Frame Size Effects...\n');
    results = run_experiment_part(x_original, t, fs, config, results, 'part_b');
end

if run_part_c
    fprintf('\nRunning Part C - Window Type Comparison...\n');
    results = run_experiment_part(x_original, t, fs, config, results, 'part_c');
end

fprintf('\nExperiment 2 completed successfully!\n');

end

% =========================================================================
% MAIN EXPERIMENT RUNNER - Eliminates duplication between parts
% =========================================================================
function results = run_experiment_part(x_original, t, fs, config, results, part_name)
% Unified function to run any experiment part with common structure

% Get part-specific configuration
part_config = get_part_config(config, part_name);

% Create log file and plots directory
[log_file, plots_dir, dir_name] = setup_output_directories(part_name, part_config);

% Run the experiment loop
signals_tested = run_experiment_loop(x_original, t, fs, part_config, log_file, plots_dir);

% Store results
results.signals = merge_structs(results.signals, signals_tested);

% Create comparison plots if multiple signals were generated
create_comparison_plots(x_original, t, fs, signals_tested, part_config, plots_dir);

% Close log file
if ~isempty(log_file) && log_file > 0
    fclose(log_file);
    fprintf('Log saved to: outputs/experiment2/%s_log.txt\n', dir_name);
end

end

% =========================================================================
% CONFIGURATION HELPER
% =========================================================================
function part_config = get_part_config(config, part_name)
% Extract part-specific configuration parameters

switch part_name
    case 'part_a'
        part_config.name = 'Hop Size Analysis';
        part_config.description = 'Testing WSOLA with different synthesis hop sizes (Hs)';
        part_config.fixed_params = struct('frame_size', config.part_a_frame_size, ...
                                        'tolerance', config.part_a_tolerance, ...
                                        'window_beta', config.part_a_window_beta, ...
                                        'alpha', config.part_a_alpha);
        part_config.test_values = config.part_a_hop_sizes;
        part_config.test_param_name = 'hop_size';
        part_config.test_param_units = 'samples';
        
    case 'part_b'
        part_config.name = 'Frame Size Effects';
        part_config.description = 'Testing different frame sizes maintaining Hs=N/2, Δ=Hs/4';
        part_config.fixed_params = struct('alpha', config.part_b_alpha, ...
                                        'window_beta', config.part_b_window_beta);
        part_config.test_values = config.part_b_frame_sizes;
        part_config.test_param_name = 'frame_size';
        part_config.test_param_units = 'samples';
        
    case 'part_c'
        part_config.name = 'Window Type Comparison';
        part_config.description = 'Testing different window types with fixed parameters';
        part_config.fixed_params = struct('alpha', config.part_c_alpha, ...
                                        'frame_size', config.part_c_frame_size, ...
                                        'syn_hop', config.part_c_syn_hop, ...
                                        'tolerance', config.part_c_tolerance);
        part_config.test_values = fieldnames(config.part_c_windows);
        part_config.test_param_name = 'window_type';
        part_config.test_param_units = '';
        part_config.windows = config.part_c_windows;
        
    otherwise
        error('Unknown part: %s', part_name);
end

end

% =========================================================================
% OUTPUT SETUP
% =========================================================================
function [log_file, plots_dir, dir_name] = setup_output_directories(part_name, part_config)
% Create output directories and log file

% Create directory structure: part_x/alpha=y.z/
alpha_val = part_config.fixed_params.alpha;
% Use more precision for alpha to handle values like 0.125
alpha_dir = sprintf('alpha=%.3f', alpha_val);

plots_dir = fullfile('outputs', 'experiment2', 'plots', part_name, alpha_dir);
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

% Use the alpha directory name for log file naming
dir_name = fullfile(part_name, alpha_dir);
log_file = fopen(fullfile('outputs', 'experiment2', sprintf('%s_log.txt', strrep(dir_name, filesep, '_'))), 'w');
if log_file > 0
    fprintf(log_file, '=== Experiment 2: %s ===\n', part_config.name);
    fprintf(log_file, 'Generated: %s\n\n', datetime("now"));
    fprintf(log_file, 'Description: %s\n\n', part_config.description);
    
    % Write fixed parameters
    fprintf(log_file, 'Fixed Parameters:\n');
    fields = fieldnames(part_config.fixed_params);
    for i = 1:length(fields)
        field = fields{i};
        value = part_config.fixed_params.(field);
        fprintf(log_file, '  %s: %s\n', field, format_value(value));
    end
    fprintf(log_file, '\n');
end

end

% =========================================================================
% MAIN EXPERIMENT LOOP
% =========================================================================
function signals_tested = run_experiment_loop(x_original, t, fs, part_config, log_file, plots_dir)
% Run the main experiment loop for any part

signals_tested = struct();

fprintf('Testing %s:\n', part_config.test_param_name);
if log_file > 0
    fprintf(log_file, '%s TESTING:\n', upper(strrep(part_config.test_param_name, '_', ' ')));
    fprintf(log_file, '%s\n\n', repmat('=', 1, length(part_config.test_param_name) + 9));
end

for i = 1:length(part_config.test_values)
    % Handle both cell arrays and numeric arrays
    if iscell(part_config.test_values)
        test_value = part_config.test_values{i};
    else
        test_value = part_config.test_values(i);
    end
    
    % Display progress
    [param_wsola, test_description, field_name] = setup_test_parameters(part_config, test_value);
    
    fprintf('  %s ... ', test_description);
    if log_file > 0
        fprintf(log_file, 'Testing %s:\n', test_description);
    end
    
    try
        % Apply WSOLA with bounds checking
        if param_wsola.alpha <= 0
            error('Alpha must be positive');
        end
        
        % Check if the hop size makes sense relative to signal length
        expected_output_length = round(length(x_original) * param_wsola.alpha);
        if expected_output_length > length(x_original) * 10
            error('Output would be too long (>10x original)');
        end
        
        y_wsola = wsolaTSM(x_original, param_wsola.alpha, param_wsola);
        
        % Check if output is reasonable
        if isempty(y_wsola) || length(y_wsola) < 100
            error('Output signal too short or empty');
        end
        
        % Store results
        signals_tested.(field_name) = y_wsola;
        
        % Generate plots
        generate_test_plots(x_original, t, y_wsola, fs, test_description, field_name, plots_dir);
        
        % Log success
        fprintf('SUCCESS\n');
        log_success(log_file, y_wsola, fs, x_original);
        
    catch ME
        % Handle failures with more detailed error info
        fprintf('FAILED: %s\n', ME.message);
        if log_file > 0
            fprintf(log_file, '  Result: FAILED - %s\n', ME.message);
            fprintf(log_file, '  Parameters: alpha=%.3f, synHop=%d, tolerance=%d, winLen=%d\n', ...
                param_wsola.alpha, param_wsola.synHop, param_wsola.tolerance, length(param_wsola.win));
            fprintf(log_file, '\n');
        end
        signals_tested.(field_name) = [];
    end
end

end

% =========================================================================
% PARAMETER SETUP
% =========================================================================
function [param_wsola, test_description, field_name] = setup_test_parameters(part_config, test_value)
% Set up WSOLA parameters for a specific test

param_wsola = struct();

switch part_config.test_param_name
    case 'hop_size'
        hop_size = test_value;
        alpha = part_config.fixed_params.alpha;
        
        param_wsola.alpha = alpha;
        param_wsola.synHop = hop_size;
        param_wsola.tolerance = part_config.fixed_params.tolerance;
        param_wsola.win = win(part_config.fixed_params.frame_size, part_config.fixed_params.window_beta);
        
        % Calculate time in milliseconds for clarity
        time_ms = hop_size / 44100 * 1000;
        test_description = sprintf('Hs=%d (%.1fms, α=%.3f)', hop_size, time_ms, alpha);
        
        field_name = sprintf('hop_%d', hop_size);
        
    case 'window_type'
        window_name = test_value;
        window_info = part_config.windows.(window_name);
        
        param_wsola.alpha = part_config.fixed_params.alpha;
        param_wsola.synHop = part_config.fixed_params.syn_hop;
        param_wsola.tolerance = part_config.fixed_params.tolerance;
        param_wsola.win = win(part_config.fixed_params.frame_size, window_info.beta);
        
        test_description = sprintf('%s window (β=%d, α=%.3f)', window_info.name, window_info.beta, part_config.fixed_params.alpha);
        field_name = sprintf('window_%s', window_name);
        
    case 'frame_size'
        frame_size = test_value;
        syn_hop = frame_size / 2;           % Hs = N/2
        tolerance = syn_hop / 4;            % Δ = Hs/4
        
        param_wsola.alpha = part_config.fixed_params.alpha;
        param_wsola.synHop = syn_hop;
        param_wsola.tolerance = tolerance;
        param_wsola.win = win(frame_size, part_config.fixed_params.window_beta);
        
        test_description = sprintf('N=%d, Hs=%d, Δ=%d, α=%.3f' , frame_size, syn_hop, tolerance, part_config.fixed_params.alpha);
        field_name = sprintf('N_%d', frame_size);
        
    otherwise
        error('Unknown test parameter: %s', part_config.test_param_name);
end

end

% =========================================================================
% PLOTTING FUNCTIONS
% =========================================================================
function generate_test_plots(x_original, t, y_wsola, fs, test_description, field_name, plots_dir)
% Generate amplitude and FFT plots for a test

% Generate amplitude vs time comparison plot
t_wsola = (0:length(y_wsola)-1) / fs;
signals_test = {x_original, y_wsola};
titles_test = {'Original (C4+G5)', test_description};

time_range = [0, 0.1];  % Show detailed time window

% Pad for comparison
if length(y_wsola) > length(x_original)
    x_padded = [x_original; zeros(length(y_wsola) - length(x_original), 1)];
    fig = plot_amplitude_vs_time({x_padded, y_wsola}, t_wsola, titles_test, time_range, false);
else
    y_padded = [y_wsola; zeros(length(x_original) - length(y_wsola), 1)];
    fig = plot_amplitude_vs_time({x_original, y_padded}, t, titles_test, time_range, false);
end

saveas(gcf, fullfile(plots_dir, sprintf('%s_comparison.png', field_name)));
close(fig);

% Generate FFT analysis
[~, ~, ~] = analyze_fft(y_wsola, fs, false);
title(sprintf('%s FFT Analysis', test_description));
saveas(gcf, fullfile(plots_dir, sprintf('%s_fft.png', field_name)));

end

function create_comparison_plots(x_original, t, fs, signals_tested, part_config, plots_dir)
% Create comparison plots of all successful results

try
    % Collect all successful results for comparison
    comparison_signals = {x_original};
    comparison_titles = {'Original'};
    
    signal_fields = fieldnames(signals_tested);
    for i = 1:length(signal_fields)
        field_name = signal_fields{i};
        if ~isempty(signals_tested.(field_name))
            comparison_signals{end+1} = signals_tested.(field_name);
            comparison_titles{end+1} = format_comparison_title(field_name, part_config);
        end
    end
    
    if length(comparison_signals) > 1
        % Find longest signal for padding
        max_length = max(cellfun(@length, comparison_signals));
        t_common = (0:max_length-1) / fs;
        
        % Pad all signals to same length
        padded_signals = cell(size(comparison_signals));
        for i = 1:length(comparison_signals)
            signal = comparison_signals{i};
            padded_signals{i} = [signal; zeros(max_length - length(signal), 1)];
        end
        
        fig = plot_amplitude_vs_time(padded_signals, t_common, comparison_titles, [0, 0.05], false);
        sgtitle(sprintf('%s Comparison', part_config.name), 'FontSize', 14, 'FontWeight', 'bold');
        
        comparison_filename = sprintf('all_%s_comparison.png', strrep(lower(part_config.test_param_name), ' ', '_'));
        saveas(gcf, fullfile(plots_dir, comparison_filename));
        close(fig);
        
        fprintf('  Comparison plot saved.\n');
    end
    
catch ME
    fprintf('  Comparison plot failed: %s\n', ME.message);
end

end

% =========================================================================
% UTILITY FUNCTIONS
% =========================================================================
function log_success(log_file, y_wsola, fs, x_original)
% Log successful test results

if log_file > 0
    fprintf(log_file, '  Result: SUCCESS\n');
    fprintf(log_file, '  Output length: %d samples (expected: ~%d)\n', ...
        length(y_wsola), length(x_original));
    
    % Quick FFT analysis for logging
    [~, ~, peaks_wsola] = analyze_fft(y_wsola, fs, false);
    if ~isempty(peaks_wsola.frequencies)
        fprintf(log_file, '  Dominant frequency: %.1f Hz\n', peaks_wsola.frequencies(1));
        fprintf(log_file, '  Peak magnitude: %.3f\n', peaks_wsola.magnitudes(1));
    end
    fprintf(log_file, '\n');
end

end

function title_str = format_comparison_title(field_name, part_config)
% Format field name into readable title for comparison plots

switch part_config.test_param_name
    case 'hop_size'
        hop_size = str2double(strrep(field_name, 'hop_', ''));
        time_ms = hop_size / 44100 * 1000;
        alpha = part_config.fixed_params.alpha;
        title_str = sprintf('Hs=%d (%.1fms)', hop_size, time_ms);
        
    case 'window_type'
        window_name = strrep(field_name, 'window_', '');
        if isfield(part_config, 'windows') && isfield(part_config.windows, window_name)
            window_info = part_config.windows.(window_name);
            alpha = part_config.fixed_params.alpha;
            title_str = sprintf('%s (β=%d, α=%.3f)', window_info.name, window_info.beta, alpha);
        else
            title_str = window_name;
        end
        
    case 'frame_size'
        frame_size = str2double(strrep(field_name, 'N_', ''));
        title_str = sprintf('N=%d (%.1fms)', frame_size, frame_size/44100*1000);
        
    otherwise
        title_str = field_name;
end

end

function value_str = format_value(value)
% Format a value for display in logs

if isnumeric(value)
    if length(value) == 1
        if isinteger(value) || value == round(value)
            value_str = sprintf('%d', value);
        else
            value_str = sprintf('%.3f', value);
        end
    else
        value_str = mat2str(value);
    end
elseif ischar(value) || isstring(value)
    value_str = char(value);
else
    value_str = class(value);
end

end

function result = merge_structs(struct1, struct2)
% Merge two structures

result = struct1;
fields = fieldnames(struct2);
for i = 1:length(fields)
    result.(fields{i}) = struct2.(fields{i});
end

end