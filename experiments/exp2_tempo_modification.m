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

% Step 3: Display alpha test values
fprintf('\nStep 3: Alpha values from configuration...\n');
fprintf('Alpha values to test: ');
fprintf('%.3f ', config.alpha_values);
fprintf('\n');

% Initialize results structure
results = struct();
results.config = config;
results.fs = fs;
results.original_signal = x_original;
results.time = t;
results.signals = struct();

% Step 4: Part A - Alpha Limits Analysis
fprintf('\nStep 4: Part A - Finding Alpha Breaking Points...\n');
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
plots_dir = fullfile('outputs', 'experiment2', 'plots','part_a');

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
        time_range = [0, 16]; 
        
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
        [~, ~, ~] = analyze_fft(x_original, fs, false);
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

fprintf('\nExperiment 2 setup completed successfully!\n');
fprintf('Ready for analysis and plotting...\n\n');

end