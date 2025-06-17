function results = exp2_tempo_modification()
% EXP2_TEMPO_MODIFICATION - Experiment 2: Tempo Modification Analysis
%
% This experiment explores WSOLA time-scale modification with focus on:
%   - Parameter boundaries (αmin = 0.125, αmax = 8)
%   - Standard parameters: N = 1024, Hs = N/2, Δ = Hs/4  
%   - Quality vs computational trade-offs
%   - Window type effects (including rectangular window limitations)
%
% Process:
%   1. Generates 2-second periodic test signal
%   2. Applies WSOLA with systematic parameter exploration
%   3. Tests extreme stretching factors (0.125x to 8x speed)
%   4. Analyzes computational performance and quality metrics
%   5. Explores window size and hop size boundaries
%
% Output:
%   results - Structure containing:
%             .signals - All processed signals at different alpha values
%             .parameters - Parameter combinations tested
%             .performance - Timing and quality metrics
%             .alpha_range - Stretching factors used
%
% Files Generated:
%   - Audio: WSOLA outputs at various stretching factors
%   - Plots: Parameter sweep results, quality vs alpha curves
%   - Performance: Computational timing analysis

fprintf('=== Experiment 2: Tempo Modification Analysis ===\n\n');

% Fixed sampling rate
fs = 44100;

% Step 1: Load configuration
fprintf('Step 1: Loading experiment configuration...\n');
config = project_config('exp2');

fprintf('  Signal duration: %.1f seconds\n', config.duration);
fprintf('  Frame size (N): %d samples\n', config.frame_size);
fprintf('  Synthesis hop (Hs): %d samples\n', config.syn_hop);
fprintf('  Tolerance (Δ): %d samples\n', config.tolerance);
fprintf('  Alpha range: %.3f to %.1f\n', config.alpha_min, config.alpha_max);

% Step 2: Generate test signal
fprintf('\nStep 2: Generating test signal...\n');
[x_original, t] = generate_periodic_signal(config.duration, config.freq1, config.freq2);

% Define output directories
output_dir = fullfile('outputs', 'experiment2');
audio_dir = fullfile(output_dir, 'audio');

fprintf('Test signal generated: %.1f seconds, %.0f Hz + %.0f Hz\n', ...
    config.duration, config.freq1, config.freq2);


% Initialize results structure
results = struct();
results.config = config;
results.fs = fs;
results.original_signal = x_original;
results.time = t;
results.signals = struct();

% Step 4: Choose which parts to run
run_part_a = false;  % Set to true to run Part A (Alpha Limits)
run_part_b = false;   % Set to true to run Part B (Window Types)
run_part_c = true;   % Set to true to run Part C (Frame Size Effects)

if run_part_a
    fprintf('\nStep 4: Running Part A - Alpha Limits Analysis...\n');
    results = run_part_a_alpha_limits(x_original, t, fs, config, results);
end

if run_part_b
    fprintf('\nStep 4: Running Part B - Window Type Comparison...\n');
    results = run_part_b_window_comparison(x_original, t, fs, config, results);
end

if run_part_c
    fprintf('\nStep 4: Running Part C - Frame Size Effects...\n');
    results = run_part_c_frame_size_effects(x_original, t, fs, config, results);
end

fprintf('\nExperiment 2 completed successfully!\n');
fprintf('Ready for analysis and plotting...\n\n');

end

function results = run_part_a_alpha_limits(x_original, t, fs, config, results)
% RUN_PART_A_ALPHA_LIMITS - Part A: Finding Alpha Breaking Points
%
% This function tests WSOLA with fixed parameters across different alpha values
% to find where the algorithm breaks or produces artifacts.

% Display alpha test values
fprintf('\nStep 3: Alpha values from configuration...\n');
fprintf('Alpha values to test: ');
fprintf('%.3f ', config.alpha_values);
fprintf('\n');
fprintf('Part A: Finding Alpha Breaking Points...\n');
fprintf('Fixed parameters: N=%d, Hs=%d, Δ=%d, Hann window\n', ...
    config.frame_size, config.syn_hop, config.tolerance);

% Create log file for Part A results
log_file = fullfile('outputs', 'experiment2', 'part_a_alpha_limits_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 2: Part A - Alpha Limits Analysis ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));
fprintf(fid, 'Fixed Parameters:\n');
fprintf(fid, '  Frame Size (N): %d samples\n', config.frame_size);
fprintf(fid, '  Synthesis Hop (Hs): %d samples\n', config.syn_hop);
fprintf(fid, '  Tolerance (Δ): %d samples\n', config.tolerance);
fprintf(fid, '  Window: Hann (beta=2)\n\n');
fprintf(fid, 'ALPHA LIMITS TESTING:\n');
fprintf(fid, '====================\n\n');

% Set up WSOLA parameters (FIXED for Part A)
param_wsola.synHop = config.syn_hop;     % Fixed: 512
param_wsola.tolerance = config.tolerance; % Fixed: 128  
param_wsola.win = win(config.frame_size, 2);  % Fixed: Hann window (beta=2)

% Create plots directory for Part A
plots_dir = fullfile('outputs', 'experiment2', 'plots', 'part_a');

fprintf('Testing alpha values for breaking points:\n');

for i = 1:length(config.alpha_values)
    alpha = config.alpha_values(i);
    
    fprintf('  Alpha = %.3f (', alpha);
    if alpha < 1
        fprintf('%.1fx faster', 1/alpha);
    elseif alpha > 1
        fprintf('%.1fx slower', alpha);
    else
        fprintf('original speed');
    end
    fprintf(') ... ');
    
    fprintf(fid, 'Testing Alpha = %.3f (', alpha);
    if alpha < 1
        fprintf(fid, '%.1fx faster', 1/alpha);
    elseif alpha > 1
        fprintf(fid, '%.1fx slower', alpha);
    else
        fprintf(fid, 'original speed');
    end
    fprintf(fid, '):\n');
    
    try
        % Apply WSOLA
        y_wsola = wsolaTSM(x_original, alpha, param_wsola);
        
        % Store results
        alpha_str = sprintf('alpha_%.3f', alpha);
        alpha_str = strrep(alpha_str, '.', 'p');  % Replace . with p for field names
        
        results.signals.(alpha_str) = y_wsola;
        
        % Generate amplitude vs time comparison plot
        t_wsola = (0:length(y_wsola)-1) / fs;  % Time vector for processed signal
        
        % Simple comparison - show both signals in a reasonable time window
        signals_alpha = {x_original, y_wsola};
        titles_alpha = {'Original (C4+G5)', sprintf('WSOLA α=%.3f', alpha)};
        
        % Show a longer time window to see tempo changes
        time_range = [0, 0.5];  % Show first 0.5 seconds
        
        % Use original time vector as reference, pad shorter signal
        if length(y_wsola) < length(x_original)
            % Processed signal is shorter (faster tempo)
            y_padded = [y_wsola; zeros(length(x_original) - length(y_wsola), 1)];
            fig = plot_amplitude_vs_time({x_original, y_padded}, t, titles_alpha, time_range, false);
        else
            % Processed signal is longer (slower tempo)  
            x_padded = [x_original; zeros(length(y_wsola) - length(x_original), 1)];
            fig = plot_amplitude_vs_time({x_padded, y_wsola}, t_wsola, titles_alpha, time_range, false);
        end
        
        saveas(gcf, fullfile(plots_dir, sprintf('part_a_alpha_%.3f_comparison.png', alpha)));
        close(fig);
        
        % Generate FFT analysis
        [~, ~, peaks_orig] = analyze_fft(x_original, fs, false);
        [~, ~, peaks_wsola] = analyze_fft(y_wsola, fs, false);
        title(sprintf('WSOLA α=%.3f FFT Analysis', alpha));
        saveas(gcf, fullfile(plots_dir, sprintf('part_a_alpha_%.3f_fft.png', alpha)));
        
        fprintf('SUCCESS\n');
        fprintf(fid, '  Result: SUCCESS\n');
        fprintf(fid, '  Output length: %d samples (expected: ~%d)\n', length(y_wsola), round(length(x_original) * alpha));
        if ~isempty(peaks_wsola.frequencies)
            fprintf(fid, '  Dominant frequency: %.1f Hz\n', peaks_wsola.frequencies(1));
        end
        fprintf(fid, '\n');
        
    catch ME
        % Handle failures at extreme parameters
        fprintf('FAILED: %s\n', ME.message);
        fprintf(fid, '  Result: FAILED - %s\n\n', ME.message);
        
        alpha_str = sprintf('alpha_%.3f', alpha);
        alpha_str = strrep(alpha_str, '.', 'p');
        
        results.signals.(alpha_str) = [];
    end
end

% Close log file
fclose(fid);
fprintf('\nPart A (Alpha Limits) completed!\n');
fprintf('Log saved to: %s\n', log_file);

end

function results = run_part_b_window_comparison(x_original, t, fs, config, results)
% RUN_PART_B_WINDOW_COMPARISON - Part B: Window Type Effects
%
% This function tests different window types with fixed WSOLA parameters
% to reveal windowing artifacts, especially rectangular window limitations.

fprintf('Part B: Window Type Comparison...\n');
fprintf('Fixed parameters: α=%.1f, N=%d, Hs=%d, Δ=%d\n', ...
    config.part_b_alpha, config.part_b_frame_size, config.part_b_syn_hop, config.part_b_tolerance);

% Create log file for Part B results
log_file = fullfile('outputs', 'experiment2', 'part_b_window_comparison_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 2: Part B - Window Type Comparison ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));
fprintf(fid, 'Fixed Parameters:\n');
fprintf(fid, '  Alpha: %.1f\n', config.part_b_alpha);
fprintf(fid, '  Frame Size (N): %d samples\n', config.part_b_frame_size);
fprintf(fid, '  Synthesis Hop (Hs): %d samples\n', config.part_b_syn_hop);
fprintf(fid, '  Tolerance (Δ): %d samples\n', config.part_b_tolerance);
fprintf(fid, '\nWINDOW TYPE TESTING:\n');
fprintf(fid, '===================\n\n');

% Create plots directory for Part B
plots_dir = fullfile('outputs', 'experiment2', 'plots', 'part_b');

% Test different window types
window_tests = {'rectangular', 'hann', 'blackman'};
fprintf('Testing window types:\n');

for i = 1:length(window_tests)
    window_name = window_tests{i};
    window_info = config.part_b_windows.(window_name);
    
    fprintf('  %s window (β=%d): %s ... ', window_info.name, window_info.beta, window_info.desc);
    fprintf(fid, 'Testing %s window (β=%d): %s\n', window_info.name, window_info.beta, window_info.desc);
    
    % Set up WSOLA parameters with specific window
    param_wsola.synHop = config.part_b_syn_hop;
    param_wsola.tolerance = config.part_b_tolerance;
    param_wsola.win = win(config.part_b_frame_size, window_info.beta);
    
    try
        % Apply WSOLA with this window type
        y_wsola = wsolaTSM(x_original, config.part_b_alpha, param_wsola);
        
        % Store results
        results.signals.(sprintf('part_b_%s', window_name)) = y_wsola;
        
        % Generate comparison plot
        t_wsola = (0:length(y_wsola)-1) / fs;
        signals_window = {x_original, y_wsola};
        titles_window = {'Original (C4+G5)', sprintf('%s Window (β=%d)', window_info.name, window_info.beta)};
        
        % Show detailed time window to see windowing artifacts
        time_range = [0, 0.1];  % Short window to see artifacts clearly
        
        % Pad for comparison
        if length(y_wsola) > length(x_original)
            x_padded = [x_original; zeros(length(y_wsola) - length(x_original), 1)];
            fig = plot_amplitude_vs_time({x_padded, y_wsola}, t_wsola, titles_window, time_range, false);
        else
            y_padded = [y_wsola; zeros(length(x_original) - length(y_wsola), 1)];
            fig = plot_amplitude_vs_time({x_original, y_padded}, t, titles_window, time_range, false);
        end
        
        saveas(gcf, fullfile(plots_dir, sprintf('part_b_%s_comparison.png', window_name)));
        close(fig);
        
        % Generate FFT analysis to see spectral artifacts
        [~, ~, peaks_wsola] = analyze_fft(y_wsola, fs, false);
        title(sprintf('%s Window FFT Analysis', window_info.name));
        saveas(gcf, fullfile(plots_dir, sprintf('part_b_%s_fft.png', window_name)));
        
        fprintf('SUCCESS\n');
        fprintf(fid, '  Result: SUCCESS\n');
        fprintf(fid, '  Output length: %d samples\n', length(y_wsola));
        if ~isempty(peaks_wsola.frequencies)
            fprintf(fid, '  Dominant frequency: %.1f Hz\n', peaks_wsola.frequencies(1));
            fprintf(fid, '  Peak magnitude: %.3f\n', peaks_wsola.magnitudes(1));
        end
        fprintf(fid, '\n');
        
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        fprintf(fid, '  Result: FAILED - %s\n\n', ME.message);
        results.signals.(sprintf('part_b_%s', window_name)) = [];
    end
end

% Create comparison plot of all window types
fprintf('  Creating window comparison plot...\n');
try
    % Collect all successful results for comparison
    comparison_signals = {x_original};
    comparison_titles = {'Original'};
    
    for i = 1:length(window_tests)
        window_name = window_tests{i};
        signal_field = sprintf('part_b_%s', window_name);
        if isfield(results.signals, signal_field) && ~isempty(results.signals.(signal_field))
            comparison_signals{end+1} = results.signals.(signal_field);
            window_info = config.part_b_windows.(window_name);
            comparison_titles{end+1} = sprintf('%s (β=%d)', window_info.name, window_info.beta);
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
        sgtitle(sprintf('Window Type Comparison (α=%.1f)', config.part_b_alpha), 'FontSize', 14, 'FontWeight', 'bold');
        saveas(gcf, fullfile(plots_dir, 'part_b_all_windows_comparison.png'));
        close(fig);
        
        fprintf('  Window comparison plot saved.\n');
    end
    
catch ME
    fprintf('  Comparison plot failed: %s\n', ME.message);
    fprintf(fid, 'Comparison plot failed: %s\n', ME.message);
end
% Close log file
fclose(fid);
fprintf('\nPart B (Window Comparison) completed!\n');
fprintf('Log saved to: %s\n', log_file);
end

function results = run_part_c_frame_size_effects(x_original, t, fs, config, results)
% RUN_PART_C_FRAME_SIZE_EFFECTS - Part C: Frame Size Effects
%
% This function tests different frame sizes with fixed WSOLA parameters
% to find optimal frame size and reveal size-related artifacts.

fprintf('Part C: Frame Size Effects...\n');
fprintf('Fixed parameters: α=%.1f, Hann window (β=%d)\n', ...
    config.part_c_alpha, config.part_c_window_beta);

% Create log file for Part C results
log_file = fullfile('outputs', 'experiment2', 'part_c_frame_size_effects_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, '=== Experiment 2: Part C - Frame Size Effects ===\n');
fprintf(fid, 'Generated: %s\n\n', datetime("now"));
fprintf(fid, 'Fixed Parameters:\n');
fprintf(fid, '  Alpha: %.1f\n', config.part_c_alpha);
fprintf(fid, '  Window: Hann (beta=%d)\n', config.part_c_window_beta);
fprintf(fid, '\nFRAME SIZE TESTING:\n');
fprintf(fid, '==================\n');
fprintf(fid, 'Maintaining relationships: Hs=N/2, Δ=Hs/4\n\n');

% Create plots directory for Part C
plots_dir = fullfile('outputs', 'experiment2', 'plots', 'part_c');

% Test different frame sizes
fprintf('Testing frame sizes (maintaining Hs=N/2, Δ=Hs/4):\n');

for i = 1:length(config.part_c_frame_sizes)
    frame_size = config.part_c_frame_sizes(i);
    syn_hop = frame_size / 2;           % Hs = N/2
    tolerance = syn_hop / 4;            % Δ = Hs/4
    
    fprintf('  N=%d, Hs=%d, Δ=%d ... ', frame_size, syn_hop, tolerance);
    fprintf(fid, 'Testing N=%d (Hs=%d, Δ=%d):\n', frame_size, syn_hop, tolerance);
    
    % Set up WSOLA parameters with specific frame size
    param_wsola.synHop = syn_hop;
    param_wsola.tolerance = tolerance;
    param_wsola.win = win(frame_size, config.part_c_window_beta);
    
    try
        % Apply WSOLA with this frame size
        y_wsola = wsolaTSM(x_original, config.part_c_alpha, param_wsola);
        
        % Store results
        results.signals.(sprintf('part_c_N_%d', frame_size)) = y_wsola;
        
        % Generate comparison plot
        t_wsola = (0:length(y_wsola)-1) / fs;
        signals_frame = {x_original, y_wsola};
        titles_frame = {'Original (C4+G5)', sprintf('N=%d samples (%.1fms)', frame_size, frame_size/fs*1000)};
        
        % Show detailed time window to see frame size effects
        time_range = [0, 0.1];  % Short window to see artifacts clearly
        
        % Pad for comparison
        if length(y_wsola) > length(x_original)
            x_padded = [x_original; zeros(length(y_wsola) - length(x_original), 1)];
            fig = plot_amplitude_vs_time({x_padded, y_wsola}, t_wsola, titles_frame, time_range, false);
        else
            y_padded = [y_wsola; zeros(length(x_original) - length(y_wsola), 1)];
            fig = plot_amplitude_vs_time({x_original, y_padded}, t, titles_frame, time_range, false);
        end
        
        saveas(gcf, fullfile(plots_dir, sprintf('part_c_N_%d_comparison.png', frame_size)));
        close(fig);
        
        % Generate FFT analysis to see spectral effects
        [~, ~, peaks_wsola] = analyze_fft(y_wsola, fs, false);
        title(sprintf('N=%d Frame Size FFT Analysis', frame_size));
        saveas(gcf, fullfile(plots_dir, sprintf('part_c_N_%d_fft.png', frame_size)));
        
        fprintf('SUCCESS\n');
        fprintf(fid, '  Result: SUCCESS\n');
        fprintf(fid, '  Output length: %d samples\n', length(y_wsola));
        fprintf(fid, '  Frame duration: %.1f ms\n', frame_size/fs*1000);
        if ~isempty(peaks_wsola.frequencies)
            fprintf(fid, '  Dominant frequency: %.1f Hz\n', peaks_wsola.frequencies(1));
            fprintf(fid, '  Peak magnitude: %.3f\n', peaks_wsola.magnitudes(1));
        end
        fprintf(fid, '\n');
        
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        fprintf(fid, '  Result: FAILED - %s\n\n', ME.message);
        results.signals.(sprintf('part_c_N_%d', frame_size)) = [];
    end
end

% Create comparison plot of all frame sizes
fprintf('  Creating frame size comparison plot...\n');
try
    % Collect all successful results for comparison
    comparison_signals = {x_original};
    comparison_titles = {'Original'};
    
    for i = 1:length(config.part_c_frame_sizes)
        frame_size = config.part_c_frame_sizes(i);
        signal_field = sprintf('part_c_N_%d', frame_size);
        if isfield(results.signals, signal_field) && ~isempty(results.signals.(signal_field))
            comparison_signals{end+1} = results.signals.(signal_field);
            comparison_titles{end+1} = sprintf('N=%d (%.1fms)', frame_size, frame_size/fs*1000);
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
        sgtitle(sprintf('Frame Size Comparison (α=%.1f)', config.part_c_alpha), 'FontSize', 14, 'FontWeight', 'bold');
        saveas(gcf, fullfile(plots_dir, 'part_c_all_frame_sizes_comparison.png'));
        close(fig);
        
        fprintf('  Frame size comparison plot saved.\n');
    end
    
catch ME
    fprintf('  Comparison plot failed: %s\n', ME.message);
    fprintf(fid, 'Comparison plot failed: %s\n', ME.message);
end

% Close log file
fclose(fid);
fprintf('\nPart C (Frame Size Effects) completed!\n');
fprintf('Log saved to: %s\n', log_file);

end



