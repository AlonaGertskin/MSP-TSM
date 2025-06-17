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
% Configuration for Experiment 2: Tempo Modification Analysis

% Get global config first
config = project_config();

% Signal generation parameters
config.duration = 2;                    % 2 seconds as specified
config.freq1 = 261;                     % C4 frequency (Hz) 
config.freq2 = 783;                     % G5 frequency (Hz)

% Standard WSOLA parameters from project overview
config.frame_size = 1024;               % N = 1024
config.syn_hop = config.frame_size / 2; % Hs = N/2 = 512
config.tolerance = config.syn_hop / 4;  % Δ = Hs/4 = 128

% Alpha range for exploration
config.alpha_min = 0.125;               % Minimum stretching factor (8x faster)
config.alpha_max = 8.0;                 % Maximum stretching factor (8x slower)

% Define test points across the range
config.alpha_values = [
    0.125,                              ... 8x faster (extreme compression)
    0.25,                               ... 4x faster  
    0.5,                                ... 2x faster
    0.75,                               ... 1.33x faster
    1.0,                                ... Original speed (reference)
    1.5,                                ... 1.5x slower
    2.0,                                ... 2x slower
    4.0,                                ... 4x slower
    8.0                                 ... 8x slower (extreme stretching)
];

% Part B: Window Type Comparison parameters
config.part_b_alpha = 2;             % Fixed alpha for window comparison
config.part_b_frame_size = 1024;       % Fixed N = 1024
config.part_b_syn_hop = 512;           % Fixed Hs = N/2
config.part_b_tolerance = 128;         % Fixed Δ = Hs/4

% Window type definitions for Part B
config.part_b_windows = struct();
config.part_b_windows.rectangular = struct('beta', 0, 'name', 'Rectangular', 'desc', 'No windowing (β=0)');
config.part_b_windows.hann = struct('beta', 2, 'name', 'Hann', 'desc', 'Standard Hann window (β=2)');
config.part_b_windows.blackman = struct('beta', 4, 'name', 'Blackman-like', 'desc', 'Higher order (β=4)');

% Part C: Frame Size Effects parameters
config.part_c_alpha = 2.0;             % Fixed alpha for frame size comparison
config.part_c_window_beta = 2;         % Fixed Hann window (β=2)

% Frame sizes to test (maintaining Hs=N/2, Δ=Hs/4 relationships)
config.part_c_frame_sizes = [256, 512, 1024, 2048, 4096];
end