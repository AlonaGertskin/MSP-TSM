function results = exp4_voice_signals()
% EXP4_VOICE_SIGNALS - Experiment 4: Voice Signal Processing
%
% Analyzes voice signals and performs gender voice conversion using TSM
% while preserving tempo. Compares OLA vs WSOLA for voice processing.

fprintf('=== Experiment 4: Voice Signal Processing ===\n\n');

% Setup
config = project_config('exp4');
plots_dir = fullfile('outputs', 'experiment4', 'plots');
audio_dir = fullfile('outputs', 'experiment4', 'audio');
if ~exist(plots_dir, 'dir'), mkdir(plots_dir); end
if ~exist(audio_dir, 'dir'), mkdir(audio_dir); end

% Load and analyze voices
fprintf('Step 1: Loading and analyzing voices...\n');
[female_voice, male_voice, fs] = load_voices(config);
[f0_female, f0_male] = analyze_voices(female_voice, male_voice, fs, config, plots_dir);

% Calculate pitch shifts
cents_f2m = 1200 * log2(f0_male / f0_female);    % Female → Male
cents_m2f = 1200 * log2(f0_female / f0_male);    % Male → Female
fprintf('  Pitch shifts: F→M = %.0f cents, M→F = %.0f cents\n', cents_f2m, cents_m2f);

% Parameter sweep
fprintf('\nStep 2: Parameter sweep...\n');
[best_params_f2m, best_params_m2f] = parameter_sweep(female_voice, male_voice, fs, ...
    f0_female, f0_male, cents_f2m, cents_m2f, config, plots_dir);

% Generate final results
fprintf('\nStep 3: Generating final results...\n');
processed_voices = generate_final_results(female_voice, male_voice, fs, ...
    f0_female, f0_male, cents_f2m, cents_m2f, best_params_f2m, best_params_m2f, plots_dir, audio_dir, config);

% Store results
results = struct();
results.fs = fs;
results.f0_female = f0_female;
results.f0_male = f0_male;
results.pitch_shifts = [cents_f2m, cents_m2f];
results.best_params_f2m = best_params_f2m;
results.best_params_m2f = best_params_m2f;
results.processed_voices = processed_voices;

fprintf('\nExperiment 4 completed!\n');

end

% =========================================================================
% VOICE LOADING
% =========================================================================
function [female_voice, male_voice, fs] = load_voices(config)

% Load female voice
female_path = fullfile('audio_inputs', config.female_voice_file);
if ~exist(female_path, 'file')
    error('Female voice file not found: %s', female_path);
end
[female_voice, fs] = audioread(female_path);

% Load male voice
male_path = fullfile('audio_inputs', config.male_voice_file);
if ~exist(male_path, 'file')
    error('Male voice file not found: %s', male_path);
end
[male_voice, fs_male] = audioread(male_path);

% Handle different sampling rates
if fs ~= fs_male
    male_voice = resample(male_voice, fs, fs_male);
end

% Convert to mono
if size(female_voice, 2) > 1, female_voice = mean(female_voice, 2); end
if size(male_voice, 2) > 1, male_voice = mean(male_voice, 2); end

% Limit duration
max_samples = round(config.analysis_duration * fs);
if length(female_voice) > max_samples, female_voice = female_voice(1:max_samples); end
if length(male_voice) > max_samples, male_voice = male_voice(1:max_samples); end

fprintf('  Loaded voices: %.1f seconds at %d Hz\n', config.analysis_duration, fs);

end

% =========================================================================
% VOICE ANALYSIS
% =========================================================================
function [f0_female, f0_male] = analyze_voices(female_voice, male_voice, fs, config, plots_dir)

% Analyze female voice
fprintf('  Analyzing female voice...\n');
[~, ~, ~] = analyze_fft(female_voice, fs, false, config.fft_freq_range);
f0_female = 191;   % Based on visual inspection of your FFT and tuner app
title('Female Voice FFT');
saveas(gcf, fullfile(plots_dir, 'female_voice_fft.png'));
close(gcf);

% Analyze male voice
fprintf('  Analyzing male voice...\n');
[~, ~, ~] = analyze_fft(male_voice, fs, false, config.fft_freq_range);
f0_male = 94;   % Based on visual inspection of your FFT and tuner app
title('Male Voice FFT');
saveas(gcf, fullfile(plots_dir, 'male_voice_fft.png'));
close(gcf);

fprintf('  F0: Female = %.1f Hz, Male = %.1f Hz\n', f0_female, f0_male);

end

% =========================================================================
% PARAMETER SWEEP
% =========================================================================
function [best_params_f2m, best_params_m2f] = parameter_sweep(female_voice, male_voice, fs, ...
    f0_female, f0_male, cents_f2m, cents_m2f, config, plots_dir)

% Initialize tracking
best_error_f2m = inf;
best_error_m2f = inf;
best_params_f2m = struct();
best_params_m2f = struct();

% Create log
log_file = fullfile('outputs', 'experiment4', 'parameter_sweep_log.txt');
fid = fopen(log_file, 'w');
fprintf(fid, 'Voice Parameter Sweep\n');
fprintf(fid, 'F→M: %.0f cents, M→F: %.0f cents\n\n', cents_f2m, cents_m2f);

test_count = 0;
total_tests = length(config.frame_sizes) * length(config.syn_hop_sizes) * length(config.tolerance_values);

for frame_size = config.frame_sizes
    for syn_hop = config.syn_hop_sizes
        for tolerance = config.tolerance_values
            
            if syn_hop >= frame_size, continue; end  % Skip invalid combinations
            
            test_count = test_count + 1;
            fprintf('  Test %d/%d: frame=%d, hop=%d, tol=%d ... ', test_count, total_tests, frame_size, syn_hop, tolerance);
            
            try
                % Test parameters
                [error_f2m, error_m2f] = test_parameters(female_voice, male_voice, fs, ...
                    f0_female, f0_male, cents_f2m, cents_m2f, frame_size, syn_hop, tolerance, config);
                
                
                % Update best F→M
                if error_f2m < best_error_f2m
                    best_error_f2m = error_f2m;
                    best_params_f2m = struct('frame_size', frame_size, 'syn_hop', syn_hop, 'tolerance', tolerance);
                    fprintf('NEW F→M! ');
                end
                
                % Update best M→F
                if error_m2f < best_error_m2f
                    best_error_m2f = error_m2f;
                    best_params_m2f = struct('frame_size', frame_size, 'syn_hop', syn_hop, 'tolerance', tolerance);
                    fprintf('NEW M→F! ');
                end
                
                fprintf('F→M:%.1f, M→F:%.1f\n', error_f2m, error_m2f);
                fprintf(fid, 'f%d_h%d_t%d: F→M=%.1f, M→F=%.1f\n', frame_size, syn_hop, tolerance, error_f2m, error_m2f);
                
            catch ME
                fprintf('FAILED\n');
                fprintf(fid, 'f%d_h%d_t%d: FAILED\n', frame_size, syn_hop, tolerance);
            end
        end
    end
end

fprintf('  Best F→M: frame=%d, hop=%d, tol=%d (%.1f Hz error)\n', ...
    best_params_f2m.frame_size, best_params_f2m.syn_hop, best_params_f2m.tolerance, best_error_f2m);
fprintf('  Best M→F: frame=%d, hop=%d, tol=%d (%.1f Hz error)\n', ...
    best_params_m2f.frame_size, best_params_m2f.syn_hop, best_params_m2f.tolerance, best_error_m2f);

fclose(fid);

end

% =========================================================================
% TEST PARAMETERS
% =========================================================================
function [error_f2m, error_m2f] = test_parameters(female_voice, male_voice, fs, ...
    f0_female, f0_male, cents_f2m, cents_m2f, frame_size, syn_hop, tolerance, config)

% Setup TSM parameters
param = struct();
param.fsAudio = fs;
param.algTSM = @wsolaTSM;
param.tolerance = tolerance;
param.synHop = syn_hop;
param.win = win(frame_size, 2);

% Apply pitch shifts
female_to_male = pitchShiftViaTSM(female_voice, cents_f2m, param);
male_to_female = pitchShiftViaTSM(male_voice, cents_m2f, param);

% Analyze results
[~, ~, peaks_f2m] = analyze_fft(female_to_male, fs, false, config.fft_freq_range);
close(gcf);
[~, ~, peaks_m2f] = analyze_fft(male_to_female, fs, false, config.fft_freq_range);
close(gcf);

% F2M error calculation with fundamental frequency detection
if ~isempty(peaks_f2m.frequencies)
    % Look for peaks near the expected fundamental (80-120 Hz for male target)
    fund_candidates = peaks_f2m.frequencies(peaks_f2m.frequencies >= 80 & peaks_f2m.frequencies <= 120);
    if ~isempty(fund_candidates)
        actual_f0_f2m = fund_candidates(1); % Use the lowest one in range
    else
        actual_f0_f2m = peaks_f2m.frequencies(1); % Fallback to highest peak
    end
    error_f2m = abs(actual_f0_f2m - 94);
else
    error_f2m = inf;
end

% M2F error calculation with fundamental frequency detection  
if ~isempty(peaks_m2f.frequencies)
    % Look for peaks near the expected fundamental (160-220 Hz for female target)
    fund_candidates = peaks_m2f.frequencies(peaks_m2f.frequencies >= 160 & peaks_m2f.frequencies <= 220);
    if ~isempty(fund_candidates)
        actual_f0_m2f = fund_candidates(1); % Use the lowest one in range
    else
        actual_f0_m2f = peaks_m2f.frequencies(1); % Fallback to highest peak
    end
    error_m2f = abs(actual_f0_m2f - 191);
else
    error_m2f = inf;
end

end

% =========================================================================
% GENERATE FINAL RESULTS
% =========================================================================
function processed_voices = generate_final_results(female_voice, male_voice, fs, ...
    f0_female, f0_male, cents_f2m, cents_m2f, best_params_f2m, best_params_m2f, plots_dir, audio_dir, config)

% Setup optimal parameters
param_f2m = struct('fsAudio', fs, 'algTSM', @wsolaTSM, 'tolerance', best_params_f2m.tolerance, ...
                   'synHop', best_params_f2m.syn_hop, 'win', win(best_params_f2m.frame_size, 2));
param_m2f = struct('fsAudio', fs, 'algTSM', @wsolaTSM, 'tolerance', best_params_m2f.tolerance, ...
                   'synHop', best_params_m2f.syn_hop, 'win', win(best_params_m2f.frame_size, 2));
param_ola = param_f2m; param_ola.tolerance = 0;  % OLA comparison

% Generate conversions
fprintf('  Applying optimal WSOLA parameters...\n');
female_to_male_wsola = pitchShiftViaTSM(female_voice, cents_f2m, param_f2m);
male_to_female_wsola = pitchShiftViaTSM(male_voice, cents_m2f, param_m2f);

fprintf('  Generating OLA comparison...\n');
female_to_male_ola = pitchShiftViaTSM(female_voice, cents_f2m, param_ola);
male_to_female_ola = pitchShiftViaTSM(male_voice, cents_m2f, param_ola);

% Create comparison plots
create_final_plots(female_voice, male_voice, female_to_male_ola, female_to_male_wsola, ...
    male_to_female_ola, male_to_female_wsola, fs, cents_f2m, cents_m2f, plots_dir);

% Save audio files
save_audio_files(female_voice, male_voice, female_to_male_ola, female_to_male_wsola, ...
    male_to_female_ola, male_to_female_wsola, fs, cents_f2m, cents_m2f, audio_dir);

% Final analysis
analyze_final_results(female_to_male_ola, female_to_male_wsola, male_to_female_ola, ...
    male_to_female_wsola, fs, f0_female, f0_male, plots_dir, config, best_params_f2m, best_params_m2f);

processed_voices = struct('f2m_ola', female_to_male_ola, 'f2m_wsola', female_to_male_wsola, ...
                         'm2f_ola', male_to_female_ola, 'm2f_wsola', male_to_female_wsola);

end

% =========================================================================
% FINAL PLOTS AND ANALYSIS
% =========================================================================
function create_final_plots(female_voice, male_voice, f2m_ola, f2m_wsola, m2f_ola, m2f_wsola, fs, cents_f2m, cents_m2f, plots_dir)

% F→M comparison
figure('Position', [100, 100, 800, 600], 'Visible', 'off');
subplot(3,1,1); plot((0:length(female_voice)-1)/fs, female_voice); title('Original Female'); xlim([0,3]);
subplot(3,1,2); plot((0:length(f2m_ola)-1)/fs, f2m_ola); title(sprintf('OLA F→M (%.0f cents)', cents_f2m)); xlim([0,3]);
subplot(3,1,3); plot((0:length(f2m_wsola)-1)/fs, f2m_wsola); title(sprintf('WSOLA F→M (%.0f cents)', cents_f2m)); xlim([0,3]);
sgtitle('Female to Male Conversion');
saveas(gcf, fullfile(plots_dir, 'final_f2m_comparison.png'));
close(gcf);

% M→F comparison  
figure('Position', [100, 100, 800, 600], 'Visible', 'off');
subplot(3,1,1); plot((0:length(male_voice)-1)/fs, male_voice); title('Original Male'); xlim([0,3]);
subplot(3,1,2); plot((0:length(m2f_ola)-1)/fs, m2f_ola); title(sprintf('OLA M→F (%.0f cents)', cents_m2f)); xlim([0,3]);
subplot(3,1,3); plot((0:length(m2f_wsola)-1)/fs, m2f_wsola); title(sprintf('WSOLA M→F (%.0f cents)', cents_m2f)); xlim([0,3]);
sgtitle('Male to Female Conversion');
saveas(gcf, fullfile(plots_dir, 'final_m2f_comparison.png'));
close(gcf);

end

function save_audio_files(female_voice, male_voice, f2m_ola, f2m_wsola, m2f_ola, m2f_wsola, fs, cents_f2m, cents_m2f, audio_dir)

audiowrite(fullfile(audio_dir, 'original_female.wav'), female_voice, fs);
audiowrite(fullfile(audio_dir, 'original_male.wav'), male_voice, fs);
audiowrite(fullfile(audio_dir, sprintf('OLA_f2m_%.0f_cents.wav', cents_f2m)), f2m_ola, fs);
audiowrite(fullfile(audio_dir, sprintf('WSOLA_f2m_%.0f_cents.wav', cents_f2m)), f2m_wsola, fs);
audiowrite(fullfile(audio_dir, sprintf('OLA_m2f_%.0f_cents.wav', cents_m2f)), m2f_ola, fs);
audiowrite(fullfile(audio_dir, sprintf('WSOLA_m2f_%.0f_cents.wav', cents_m2f)), m2f_wsola, fs);

end

function analyze_final_results(f2m_ola, f2m_wsola, m2f_ola, m2f_wsola, fs, f0_female, f0_male, plots_dir, config, best_params_f2m, best_params_m2f)
fprintf('  Final FFT analysis...\n');

% F→M analysis
[~, ~, peaks_f2m_ola] = analyze_fft(f2m_ola, fs, false, config.fft_freq_range);
title(sprintf('Final F→M OLA FFT (frame size=%d,hop size=%d)', best_params_f2m.frame_size, best_params_f2m.syn_hop));
saveas(gcf, fullfile(plots_dir, 'final_f2m_ola_fft.png')); close(gcf);

[~, ~, peaks_f2m_wsola] = analyze_fft(f2m_wsola, fs, false, config.fft_freq_range);
title(sprintf('Final F→M WSOLA FFT (frame size=%d,hop size=%d,tolerance=%d)', best_params_f2m.frame_size, best_params_f2m.syn_hop, best_params_f2m.tolerance));
saveas(gcf, fullfile(plots_dir, 'final_f2m_wsola_fft.png')); close(gcf);

% M→F analysis
[~, ~, peaks_m2f_ola] = analyze_fft(m2f_ola, fs, false, config.fft_freq_range);
title(sprintf('Final M→F OLA FFT (frame size=%d,hop size=%d)', best_params_m2f.frame_size, best_params_m2f.syn_hop));
saveas(gcf, fullfile(plots_dir, 'final_m2f_ola_fft.png')); close(gcf);

[~, ~, peaks_m2f_wsola] = analyze_fft(m2f_wsola, fs, false, config.fft_freq_range);
title(sprintf('Final M→F WSOLA FFT (frame size=%d,hop size=%d,tolerance=%d)', best_params_m2f.frame_size, best_params_m2f.syn_hop, best_params_m2f.tolerance));
saveas(gcf, fullfile(plots_dir, 'final_m2f_wsola_fft.png')); close(gcf);

% Print results
if ~isempty(peaks_f2m_ola.frequencies) && ~isempty(peaks_f2m_wsola.frequencies)
    error_ola_f2m = abs(peaks_f2m_ola.frequencies(1) - f0_male);
    error_wsola_f2m = abs(peaks_f2m_wsola.frequencies(1) - f0_male);
    fprintf('  F→M: OLA error = %.1f Hz, WSOLA error = %.1f Hz\n', error_ola_f2m, error_wsola_f2m);
end

if ~isempty(peaks_m2f_ola.frequencies) && ~isempty(peaks_m2f_wsola.frequencies)
    error_ola_m2f = abs(peaks_m2f_ola.frequencies(1) - f0_female);
    error_wsola_m2f = abs(peaks_m2f_wsola.frequencies(1) - f0_female);
    fprintf('  M→F: OLA error = %.1f Hz, WSOLA error = %.1f Hz\n', error_ola_m2f, error_wsola_m2f);
end

end