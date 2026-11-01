function fig1_paper(node_count, time, fs, outputs, E0, S0, F0, varargin)
% FIG1_PAPER  Paper-quality panel figure: time series + PSD for N parameter settings.
%
%   fig1_paper(node_count, time, fs, outputs, E0, S0, F0, 'Name', value, ...)
%
% 输入
%   node_count  列数，每个节点一列
%   time        仿真总时长（s）
%   fs          采样率（Hz）
%   outputs     node_count × n 的 EEG 矩阵（每行一个节点）
%   E0, S0, F0  三个增益；S0 若为向量则逐列标注在列标题上
%
% 可选参数（Name-Value）
%   't_skip'    跳过开头多少秒的瞬态，默认 1
%   't_view'    显示多长（s），默认 5
%   'f_max'     PSD 频率上限（Hz），默认 40
%   'SavePath'  导出路径；省略则只在屏幕上显示
%
% 输出
%   无返回值；有 SavePath 时导出矢量 PDF
%
% 做什么
%   上排：时域波形；下排：归一化 PSD，底纹标出 δ/θ/α/β/γ 频段；
%   每列顶部标题写该列用的参数值，子图编号 (a1)(a2)... 便于正文引用。
%
% 调用关系
%   上游：run_fig1.m（论文 Fig. 3）
%
% 备注
%   列数由 node_count 决定，改 run_fig1.m 里的 B 向量会自动改列数。
p = inputParser;
addParameter(p,'t_skip',   1,  @isscalar);
addParameter(p,'t_view',   5,  @isscalar);
addParameter(p,'f_max',   40,  @isscalar);
addParameter(p,'SavePath',[],  @ischar);
parse(p,varargin{:});

t_skip   = p.Results.t_skip;
t_view   = p.Results.t_view;
f_max    = p.Results.f_max;
SavePath = p.Results.SavePath;

%% ===== 全局绘图样式 =====
set(groot, ...
    'DefaultAxesFontName',     'Arial', ...
    'DefaultTextFontName',     'Arial', ...
    'DefaultAxesFontSize',     9, ...
    'DefaultLineLineWidth',    1.2);

tvec = (0:size(outputs,2)-1)/fs;

%% ===== 跳过初始瞬态 =====
t_start = t_skip;
t_end   = t_skip + t_view;
t_idx   = tvec >= t_start & tvec <= t_end;
tvec_plot = tvec(t_idx);

%% ===== Layout =====
fig = figure('Color','w','Units','centimeters');
fig.Position(3:4) = [node_count*3.2 + 2.5, 7];  % 加宽以容纳左侧编号

t = tiledlayout(fig, 2, node_count, ...
    'TileSpacing', 'tight', ...
    'Padding',     'compact');

%% ===== 节律频段 =====
bands = { ...
    'δ', 0.5, 4,   [0.00 0.45 0.74];
    'θ', 4,   8,   [0.47 0.67 0.19];
    'α', 8,  13,   [0.93 0.69 0.13];
    'β',13,  30,   [0.89 0.10 0.11];
    'γ',30,  f_max,[0.50 0.00 0.50]};

%% ===== 主循环 =====
for i = 1:node_count
    %% ---- 第一行：时域 ----
    ax1 = nexttile(t, i);
    plot(ax1, tvec_plot, outputs(i,t_idx), 'k', 'LineWidth', 1)
    
    xlabel(ax1, 'Time (s)')
    
    % 仅第一列显示y轴标签和刻度，但所有子图都显示边框
    if i == 1
        ylabel(ax1, 'Amplitude')
        ax1.YAxis.Visible = 'on';
    else
        % 关键：关闭刻度标签但保留轴线
        ax1.YAxis.TickLabels = {};
        ax1.YAxis.TickLength = [0 0];  % 隐藏刻度线
        ax1.YAxis.Color = 'k';         % 确保轴线可见
    end
    
    % 标题仅显示G_s值
    title(ax1, sprintf('$G_s = %.2f$', S0(i)), ...
        'FontSize', 9, 'FontWeight', 'normal', ...
        'Interpreter', 'latex')
    
    % 添加子图编号 (a), (b), (c)... 放在图外左上角
    pos1 = ax1.Position;  % [left bottom width height] 归一化坐标
    
    % 编号位置：子图左边界左侧，上边界上方
    annotation('textbox', ...
        [pos1(1), pos1(2)+pos1(4)+0.020, 0.02, 0.02], ...
        'String', ['(', char('a'+i-1), ')'], ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none');
    
    grid(ax1,'on')
    xlim(ax1, [t_start t_end])
    ylim(ax1, [-5 15])
    box(ax1,'on')                    % 确保完整边框
    set(ax1,'Layer','top')
    set(ax1,'TickDir','out')         % 刻度朝外，更美观
    
    %% ---- 第二行：PSD ----
    ax2 = nexttile(t, node_count + i);
    
    x = outputs(i,t_idx) - mean(outputs(i,t_idx));
    win = min(1024, length(x));
    noverlap = floor(win/2);
    nfft = max(1024, 2^nextpow2(win));
    [pxx, f] = pwelch(x, hanning(win), noverlap, nfft, fs);
    
    idx = f >= 1 & f <= f_max;
    pxx = pxx(idx); f = f(idx);
    pxx = pxx / max(pxx);
    
    hold(ax2,'on')
    for b = 1:size(bands,1)
        bs = bands{b,2}; be = min(bands{b,3}, f_max); c = bands{b,4};
        if be > bs
            h = fill(ax2, [bs bs be be], [0 1.2 1.2 0], c, ...
                'FaceAlpha',0.12,'EdgeColor','none');
            uistack(h,'bottom')
            text(ax2, mean([bs be]), 1.1, bands{b,1}, ...
                'FontSize',7,'FontWeight','bold', ...
                'HorizontalAlignment','center')
        end
    end
    plot(ax2, f, pxx, 'k', 'LineWidth', 1)
    hold(ax2,'off')
    
    xlabel(ax2, 'Frequency (Hz)')
    
    % 仅第一列显示y轴标签和刻度，但所有子图都显示边框
    if i == 1
        ylabel(ax2, 'Normalized power')
        ax2.YAxis.Visible = 'on';
    else
        % 关键：关闭刻度标签但保留轴线
        ax2.YAxis.TickLabels = {};
        ax2.YAxis.TickLength = [0 0];  % 隐藏刻度线
        ax2.YAxis.Color = 'k';         % 确保轴线可见
    end
    
    grid(ax2,'on')
    xlim(ax2, [1 f_max])
    ylim(ax2, [0 1.2])
    box(ax2,'on')                    % 确保完整边框
    set(ax2,'Layer','top')
    set(ax2,'TickDir','out')         % 刻度朝外，更美观
end

%% ===== 导出 =====
% 注意（R2024b 已知坑）：如果图窗自始至终没有显示过，exportgraphics 偶尔会把
% 各个 axes 的交互工具栏（放大镜/手掌/房子那一排图标）当成图元渲染进 PDF，
% 而且是否出现是随机的——同一份数据两次导出可能一次带图标、一次不带。
% 这里在导出前显式关掉 fig 和所有 axes 的交互元素，保证输出稳定。
try
    fig.ToolBar      = 'none';
    fig.MenuBar      = 'none';
    ax_all = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax_all)
        try, ax_all(k).Toolbar.Visible = 'off'; catch, end
        try, disableDefaultInteractivity(ax_all(k)); catch, end
    end
catch
end

if ~isempty(SavePath)
    exportgraphics(fig, SavePath, ...
        'ContentType','vector', ...
        'BackgroundColor','none', ...
        'Resolution',600)
end

end