function fig = plot_amplitude_vs_time(signals, time, titles, time_range, show_plot)
% PLOT_AMPLITUDE_VS_TIME - Plot amplitude vs time for multiple signals
%
% Input:
%   signals    - Cell array of signal vectors
%   time       - Time vector
%   titles     - Cell array of signal titles
%   time_range - [start, end] time range to plot (optional, default: full range)
%   show_plot  - Boolean to show plot (optional, default: false)
%
% Output:
%   fig        - Figure handle

% Set defaults
if nargin < 4 || isempty(time_range)
    time_range = [time(1), time(end)];
end
if nargin < 5
    show_plot = false;
end

% Find indices for time range
start_idx = find(time >= time_range(1), 1, 'first');
end_idx = find(time <= time_range(2), 1, 'last');

% Create subplot layout
num_signals = length(signals);
if show_plot
    fig = figure('Position', [100, 100, 800, 150*num_signals]);
else
    fig = figure('Position', [100, 100, 800, 150*num_signals], 'Visible', 'off');
end

for i = 1:num_signals
    subplot(num_signals, 1, i);
    
    % Plot signal in specified time range
    plot(time(start_idx:end_idx), signals{i}(start_idx:end_idx), 'b-', 'LineWidth', 1);
    
    % Format plot
    grid on;
    title(titles{i}, 'FontSize', 12);
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim(time_range);
    
    % Set y-limits for better visualization
    y_max = max(abs(signals{i}(start_idx:end_idx)));
    if y_max > 0
        ylim([-y_max*1.1, y_max*1.1]);
    end
end

% Add overall title
sgtitle('Amplitude vs Time Comparison', 'FontSize', 14, 'FontWeight', 'bold');

end