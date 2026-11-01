% ==============================
% plaintext: 长时程分类结果可视化（修复图例对齐问题）
% 上子图：原始波形，下子图：类型色带，右侧独立图例不挤占主图空间
% ==============================
function plot_timeseries_result(signal, Fs, t_type, type_seq, params)
signal = signal(:);
t_signal = (0:length(signal)-1)/Fs;
dt = mean(diff(t_type)); % 步长，用于色块宽度

%% ----------------------
% 创建图形，手动设置子图位置（像素级对齐）
% ----------------------
fig = figure('Position', [100 100 1600 520], 'Color', 'white');

% 上子图：原始波形
ax1 = axes('Parent', fig, 'Position', [0.06 0.56 0.88 0.38]);
plot(ax1, t_signal, signal, 'k', 'LineWidth', 0.8);
xlabel(ax1, 'Time (s)', 'FontSize', 11);
ylabel(ax1, 'Amplitude (mV)', 'FontSize', 11);
title(ax1, 'Raw EEG Signal', 'FontSize', 13, 'FontWeight', 'bold');
grid(ax1, 'on'); box(ax1, 'on');
xlim(ax1, [t_signal(1), t_signal(end)]);
set(ax1, 'TickDir', 'out');

% 下子图：类型色带（和上方等宽，像素级对齐）
ax2 = axes('Parent', fig, 'Position', [0.06 0.08 0.88 0.38]);
hold(ax2, 'on');
for k = 1:length(t_type)
    t_start = t_type(k) - dt/2;
    t_end = t_type(k) + dt/2;
    c = params.type_colors{type_seq(k)};
    patch(ax2, [t_start t_end t_end t_start], [0 0 1 1], c, ...
        'EdgeColor', 'none', 'FaceAlpha', 1);
end
xlabel(ax2, 'Time (s)', 'FontSize', 11);
ylabel(ax2, 'Type', 'FontSize', 11);
title(ax2, 'Classified Signal Type (Color Block Width = Step Size)', 'FontSize', 13, 'FontWeight', 'bold');
ylim(ax2, [0 1]);
set(ax2, 'YTick', [], 'TickDir', 'out');
grid(ax2, 'on'); box(ax2, 'on');
xlim(ax2, [t_signal(1), t_signal(end)]);

%% ----------------------
% 右侧独立图例（完全不干扰主图对齐）
% ----------------------
leg_ax = axes('Parent', fig, 'Position', [0.95 0.15 0.04 0.7]);
hold(leg_ax, 'on');
for i = 1:6
    patch(leg_ax, [0 1 1 0], [i-1 i-1 i i]/6.5, params.type_colors{i}, ...
        'EdgeColor', 'k', 'LineWidth', 0.5);
    text(leg_ax, 1.3, (i-0.5)/6.5, params.type_names{i}, 'FontSize', 8, ...
        'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end
set(leg_ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
xlim(leg_ax, [0 3]); ylim(leg_ax, [0 1]);
set(leg_ax, 'TickDir', 'out');

%% ----------------------
% 对齐上下子图x轴
% ----------------------
linkaxes([ax1, ax2], 'x');
end