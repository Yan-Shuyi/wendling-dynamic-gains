function fig = fig3_paper(dy_fs, dy_fs_srate, time_points, period_labels, out_path)
% FIG3_PAPER  Boxplots of E/S/F gains across experiments and seizure periods.
%
%   fig = fig3_paper(dy_fs, dy_fs_srate, time_points, period_labels, out_path)
%
% 输入
%   dy_fs          3*N × n 的增益轨迹（N 次实验，每实验 3 行：G_e, G_s, G_f）
%   dy_fs_srate    增益轨迹采样率（Hz）
%   time_points    四个取样时刻（s），默认 [10 30 50 70]
%   period_labels  四个时期的标签，默认 interICTAL/preONSET/ONSET/ICTAL
%   out_path       导出路径；为空则只返回图句柄不导出
%
% 输出
%   fig    figure 句柄
%
% 做什么
%   三个子图分别画 G_e（EXC）、G_s（SDI）、G_f（FSI）；
%   每个子图里四个箱子对应四个时期，黑白发表风格（参考原文 Figure 7b）。
%   用分布而非单条轨迹说明"四个时期的增益是可分的"。
%
% 调用关系
%   上游：run_fig3.m（论文 Fig. 5）
%
% 备注
%   N 由 size(dy_fs,1)/3 推出，即"实验次数"；因为 sc_matrix 全零，
%   N 次实验之间没有任何耦合，等价于 N 次独立重复。
    if nargin < 4 || isempty(period_labels)
        period_labels = {'interICTAL', 'preONSET', 'ONSET', 'ICTAL'};
    end
    if nargin < 3 || isempty(time_points)
        time_points = [10, 30, 50, 70];
    end
    if nargin < 5 || isempty(out_path)
        out_path = '';
    end

    num_periods = length(time_points);
    N = size(dy_fs, 1) / 3;

    E_mat = zeros(N, num_periods);
    S_mat = zeros(N, num_periods);
    F_mat = zeros(N, num_periods);

    for j = 1:N
        E = dy_fs(3*j-2, :);
        S = dy_fs(3*j-1, :);
        F = dy_fs(3*j, :);

        t_dy = (0:length(E)-1) / dy_fs_srate;

        for p = 1:num_periods
            [~, idx] = min(abs(t_dy - time_points(p)));
            E_mat(j, p) = E(idx);
            S_mat(j, p) = S(idx);
            F_mat(j, p) = F(idx);
        end
    end
    % ========== 数据部分结束 ==========

    % ========== 出图部分：手动定位 ==========
    set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
             'DefaultTextFontName', 'Times New Roman');

    fig = figure('Name', 'ESF Boxplot (Multi-run)', 'Color', 'w', ...
                 'Position', [100, 100, 1050, 380]);
    fig.PaperPositionMode = 'auto';

    % ---- 手动指定三个子图的位置 [left, bottom, width, height] ----
    % 归一化坐标 (0~1)，相对于整个 figure
    % 布局：左边留 0.07，右边留 0.03，三个子图各占 0.28，间距 0.04
    % 底部留 0.14（给 X 标签），顶部留 0.12（给 (a)(b)(c) 标签）

    ax1 = axes('Position', [0.07, 0.14, 0.28, 0.74]);
    ax2 = axes('Position', [0.39, 0.14, 0.28, 0.74]);
    ax3 = axes('Position', [0.71, 0.14, 0.28, 0.74]);

    %% ---- 子图1: EXC ----
    h1 = boxplot(ax1, E_mat, 'Labels', period_labels, ...
                 'Colors', 'k', 'Symbol', 'k+', 'Widths', 0.6);
    set(h1, {'Color'}, {'k'});
    ylabel(ax1, 'EXC', 'FontSize', 11);
    set(ax1, 'FontSize', 10, 'LineWidth', 1);
    grid(ax1, 'off');

    %% ---- 子图2: SDI ----
    h2 = boxplot(ax2, S_mat, 'Labels', period_labels, ...
                 'Colors', 'k', 'Symbol', 'k+', 'Widths', 0.6);
    set(h2, {'Color'}, {'k'});
    ylabel(ax2, 'SDI', 'FontSize', 11);
    set(ax2, 'FontSize', 10, 'LineWidth', 1);
    grid(ax2, 'off');

    %% ---- 子图3: FSI ----
    h3 = boxplot(ax3, F_mat, 'Labels', period_labels, ...
                 'Colors', 'k', 'Symbol', 'k+', 'Widths', 0.6);
    set(h3, {'Color'}, {'k'});
    ylabel(ax3, 'FSI', 'FontSize', 11);
    set(ax3, 'FontSize', 10, 'LineWidth', 1);
    grid(ax3, 'off');

    %% ---- 用 annotation 放置 (a)(b)(c) ----
    annotation(fig, 'textbox', [0.05, 0.90, 0.04, 0.04], ...
        'String', '(a)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    annotation(fig, 'textbox', [0.37, 0.90, 0.04, 0.04], ...
        'String', '(b)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    annotation(fig, 'textbox', [0.69, 0.90, 0.04, 0.04], ...
        'String', '(c)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    drawnow;

    %% ---- 导出 ----
    if ~isempty(out_path)
        exportgraphics(fig, out_path, ...
            'Resolution', 300, ...
            'ContentType', 'vector');
    end

end