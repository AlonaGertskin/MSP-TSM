function config = project_config(experiment)
% PROJECT_CONFIG - Configuration parameters for TSM project
%
% Usage:
%   config = project_config()           % Returns global config
%   config = project_config('exp1')     % Returns experiment 1 config
%   config = project_config('exp2')     % Returns experiment 2 config

if nargin == 0
    % Global parameters shared across all experiments
    config.fs = 44100;                % Default sampling rate (Hz)
    config.output_dir = 'outputs';    % Output directory
else
    switch experiment
        case 'exp1'
            config = config_exp1();
        case 'exp1_sweep'
            config = config_exp1_sweep();
        case 'exp2'
            config = config_exp2();
        otherwise
            error('Unknown experiment: %s', experiment);
    end
end

end

function config = config_exp1()
% Configuration for Experiment 1: Pitch Modification with Periodic Signals

% Get global config first
config = project_config();

% Add experiment-specific parameters
config.duration = 1;                  % Signal duration (seconds)
config.freq1 = 261;                   % C4 frequency (Hz)
config.freq2 = 783;                   % G5 frequency (Hz) - 3rd harmonic
config.amp1 = 1.0;                    % Primary amplitude
config.amp2 = 0.4;                    % Secondary amplitude
config.target_g4 = 392;               % G4 (Dominant)
config.target_f3 = 175;               % F3 (Sub-Dominant)

end

function config = config_exp1_sweep()
% Configuration for Experiment 1 parameter sweep

% Parameter sweep values for OLA
config.ola_synHop_values = [512,1024];
config.ola_winLen_values = [256,1024];
config.ola_winType_values = [2];

%Parameter sweep values for WSOLA
config.wsola_tolerance_values = [512, 1024, 2048, 4096];

end

function config = config_exp2()
% Configuration for Experiment 2: Tempo Modification Analysis (Updated Structure)

% Get global config first
config = project_config();

% Signal generation parameters (shared across all parts)
config.duration = 2;                    % 2 seconds as specified
config.freq1 = 261;                     % C4 frequency (Hz) 
config.freq2 = 783;                     % G5 frequency (Hz)

% Standard base parameters for all parts
config.frame_size = 1024;               % Standard frame size for Parts A and C
config.tolerance = 128;                 % Standard tolerance (frame_size/8)

% ==========================================================================
% PART A: Hop Size Analysis
% ==========================================================================
config.part_a_alpha = 0.125;             % Fixed alpha for hop size testing
config.part_a_frame_size = 1024;       % Fixed frame size
config.part_a_tolerance = 128;         % Fixed tolerance
config.part_a_window_beta = 2;         % Fixed Hann window

% Hop size values to test (synthesis hop sizes in samples)
config.part_a_hop_sizes = [128, 256, 512, 1024, 2048];

% ==========================================================================
% PART B: Frame Size Effects  
% ==========================================================================
config.part_b_alpha = 0.125;             % Fixed alpha for frame size comparison
config.part_b_window_beta = 2;         % Fixed Hann window (β=2)

% Frame sizes to test (maintaining Hs=N/2, Δ=Hs/4 relationships)
config.part_b_frame_sizes = [256, 512, 1024, 2048, 4096];

% ==========================================================================
% PART C: Window Type Comparison
% ==========================================================================
config.part_c_alpha = 0.125;             % Fixed alpha for window comparison
config.part_c_frame_size = 1024;       % Fixed N = 1024
config.part_c_syn_hop = 512;           % Fixed Hs = N/2
config.part_c_tolerance = 128;         % Fixed Δ = Hs/4

% Window type definitions for Part C
config.part_c_windows = struct();
config.part_c_windows.rectangular = struct('beta', 0, 'name', 'Rectangular', 'desc', 'No windowing (β=0)');
config.part_c_windows.hann = struct('beta', 2, 'name', 'Hann', 'desc', 'Standard Hann window (β=2)');
config.part_c_windows.sin_4_window = struct('beta', 4, 'name', 'sin^4 window', 'desc', 'Higher order (β=4)');

end