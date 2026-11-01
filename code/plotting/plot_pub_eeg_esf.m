% ==============================
% plaintext: 论文级EEG-ESF联合可视化函数（最终稳定版）
% 修复所有已知Bug：
% 1. patch LineWidth必须>0（MATLAB硬性要求）
% 2. 时间轴直接使用原始向量，不做任何清洗和转置
% 3. 类型区间检测用diff+find，逻辑最简洁
% ==============================
function plot_pub_eeg_esf(output, type_full, Ge, Gs, Gf, Fs, params, dy_fs)
% PLOT_PUB_EEG_ESF  Three-layer paper figure: EEG / state-type bar / gain traces.
%
%   plot_pub_eeg_esf(output, type_full, Ge, Gs, Gf, Fs, params, dy_fs)
%
% 输入
%   output     整段 EEG（1 × N，mV）
%   type_full  逐采样点的类别序列（classify_eeg_timeseries 的第三个输出）
%   Ge,Gs,Gf   三条增益轨迹（1 × n_g）
%   Fs         EEG 采样率（Hz）
%   params     classify_params()，用来取六类的颜色/名称
%   dy_fs      增益轨迹采样率（Hz）
%
% 输出
%   无返回值；只画到当前 figure，导出由调用方用 exportgraphics 完成
%
% 做什么
%   上：EEG 波形（线宽 0.2，200 s 压缩后仍能看出振幅包络）
%   中：六类状态条，每个色块上标类型编号——这条带子就是论文 Fig. 2 那棵树的输出
%   下：G_e / G_s / G_f 三条轨迹，用双 y 轴与上层对齐时间
%
% 调用关系
%   上游：run_fig2.m（论文 Fig. 4）
%
% 备注
%   1. exportgraphics 的紧裁切会切掉 y 轴标签字头（实测标签落在页面外
%      xMin = -4.15 pt）。这里在标签左侧约 2.6% 轴宽处放一个白色竖排字符，
%      把裁切边界推到标签之外——这是试过多种方案后唯一有效且不降级字体的做法
%      （print -dpdf 会把字体降级成未嵌入的 base-14，期刊不接受）。
%   2. patch 的 LineWidth 必须 > 0，这是 MATLAB 的硬性要求。
output = output(:); type_full = type_full(:);
Ge = Ge(:); Gs = Gs(:); Gf = Gf(:);
N = length(output);

% ESF重采样
if nargin < 8 || isempty(dy_fs), dy_fs = Fs; end
t_out = (0:N-1)/Fs;
if dy_fs ~= Fs
    t_esf = (0:length(Ge)-1)/dy_fs;
    Ge = interp1(t_esf, Ge, t_out, 'linear', 'extrap');
    Gs = interp1(t_esf, Gs, t_out, 'linear', 'extrap');
    Gf = interp1(t_esf, Gf, t_out, 'linear', 'extrap');
end

% ----------------------
% 样式定义（LineWidth全部>=0.5，满足MATLAB要求）（方案A：类型1/2差异化，黑白打印友好）
% ----------------------
if ~isfield(params, 'type_styles')
    % Type 1: 背景 → 白底+极细黑实线边框（解决“空洞感”）
    params.type_styles(1) = struct('fc',[1 1 1], 'ec',[0 0 0], 'ls','-', 'lw',0.5);
    % Type 2: 零星棘波 → 极浅灰底+细点线边框（和Type1形成区分）
    params.type_styles(2) = struct('fc',[0.92 0.92 0.92], 'ec',[0 0 0], 'ls',':', 'lw',0.6);
    % Type 3: 持续棘波 → 浅灰底+黑细实线（保持不变，和Type2形成梯度）
    params.type_styles(3) = struct('fc',[0.9 0.9 0.9], 'ec',[0 0 0], 'ls','-', 'lw',0.6);
    % Type 4: 慢节律 → 中灰底+黑实线（保持不变）
    params.type_styles(4) = struct('fc',[0.7 0.7 0.7], 'ec',[0 0 0], 'ls','-', 'lw',1.0);
    % Type 5: 快活动 → 深灰底+黑点划线（保持不变）
    params.type_styles(5) = struct('fc',[0.4 0.4 0.4], 'ec',[0 0 0], 'ls','-.','lw',1.0);
    % Type 6: 高幅慢波 → 纯黑（保持不变，最醒目）
    params.type_styles(6) = struct('fc',[0 0 0], 'ec','none', 'ls','-', 'lw',0.5);
end
if ~isfield(params,'mid_layer_height'), params.mid_layer_height = 0.05; end
if ~isfield(params,'output_ylim'), params.output_ylim = [-7, 25]; end
if ~isfield(params,'esf_ylim'), params.esf_ylim = [2.5, 30.5]; end

%% ----------------------
% 1. 创建画布（固定像素级间隙，确保中层/下层有呼吸空间）
% ----------------------
fig = figure('Color','white','Position',[100 100 800 400]);
lm = 0.115; rm = 0.03; w = 1 - lm - rm;   % lm 需容纳 11pt 粗体 y 轴标签，0.09 会裁字

% 固定间隙（归一化单位，约占总高的1.5%）
gap = 0.015;

% 中层高度：原为 0.03，打印后类型条几乎不可辨；改为 0.085 以容纳类型编号
mid_h = 0.085;

% 上下两层均分剩余高度
remaining = 1 - 0.03 - 0.09 - mid_h - 2*gap; % top + bottom + mid + 2 gaps
top_h = remaining * 0.4;  % 上层40%
bot_h = remaining * 0.6;  % 下层60%（ESF需要更多空间）

% 从上到下排列
b1 = 1 - 0.03 - top_h;        % 上层底部
h1 = top_h;

mb = b1 - gap - mid_h;        % 中层底部
mh = mid_h;

b3 = mb - gap - bot_h;        % 下层底部
h3 = bot_h;

% 创建三个Axes
ax1 = axes('Parent',fig,'Position',[lm b1 w h1]);  % 上层EEG
ax2 = axes('Parent',fig,'Position',[lm mb w mh]);  % 中层Type Bar
ax3 = axes('Parent',fig,'Position',[lm b3 w h3]);  % 下层ESF

linkaxes([ax1,ax2,ax3],'x');

%% ----------------------
% 2. 上层：EEG
% ----------------------
% 200 s × 256 Hz 共 5 万余点画在约 12 cm 宽度内；线宽必须很细，
% 否则波形会糊成一整块黑色，看不到振幅包络的起伏。
plot(ax1, t_out, output, 'k-', 'LineWidth', 0.2, 'Color', [0.15 0.15 0.15]);
ylabel(ax1,'Amplitude (mV)','FontName','Times New Roman','FontSize',11,'FontWeight','bold');
set(ax1,'XTickLabel',[],'YGrid','on','Box','on','YLim',params.output_ylim);
set(ax1,'FontName','Times New Roman','FontSize',10);

%% ----------------------
% 3. 中层：Type Bar（标准游程编码，无任何数据清洗）
% ----------------------
% 找变化点：diff(type_full)非零的位置
d = diff(type_full);
cp = find(d ~= 0); % change points

% 构建区间
seg_start = [1; cp+1];
seg_end = [cp; length(type_full)];

fprintf('共检测到 %d 个类型区间，开始绘制...\n', length(seg_start));

for k = 1:length(seg_start)
    s = seg_start(k); e = seg_end(k);
    ts = t_out(s); te = t_out(e);
    ct = type_full(s);
    
    % 跳过无效
    if ts >= te, continue; end
    if ct < 1 || ct > 6, continue; end
    
    st = params.type_styles(ct);
    
    % 【关键】LineWidth必须>0，MATLAB硬性要求
    lw = max(st.lw, 0.01);
    
    % 当EdgeColor为'none'时，边框不显示，LineWidth无所谓
    % 当EdgeColor为黑色时，LineWidth控制边框粗细
    patch(ax2, [ts te te ts], [0 0 1 1], st.fc, ...
        'EdgeColor', st.ec, ...
        'LineStyle', st.ls, ...
        'LineWidth', lw, ...
        'FaceAlpha', 1);

    % 在足够宽的区间内标注类型编号（配合论文中 Table 的六类定义）
    if (te - ts) > 6 && ct <= 3        % 浅底区间用黑字
        text(ax2, (ts+te)/2, 0.5, sprintf('%d', ct), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontName','Times New Roman','FontSize',8,'Color','k');
    elseif (te - ts) > 6              % 深底区间用白字
        text(ax2, (ts+te)/2, 0.5, sprintf('%d', ct), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontName','Times New Roman','FontSize',8,'Color','w');
    end
end

set(ax2,'XTick',[],'YTick',[],'XColor','none','YColor','none','Box','off');
set(ax2,'YLim',[0 1],'XLim',[t_out(1), t_out(end)]);

%% ----------------------
% 4. 下层：ESF
% ----------------------
plot(ax3, t_out, Ge, 'k-',  'LineWidth',1.5,'DisplayName','$G_e$'); hold(ax3,'on');
plot(ax3, t_out, Gs, 'k-.', 'LineWidth',1.0,'DisplayName','$G_s$');
plot(ax3, t_out, Gf, 'k--', 'LineWidth',1.0,'DisplayName','$G_f$'); hold(ax3,'off');

xlabel(ax3,'Time (s)','FontName','Times New Roman','FontSize',11,'FontWeight','bold');
ylabel(ax3,'Synaptic gain (mV)','FontName','Times New Roman','FontSize',11,'FontWeight','bold');
legend(ax3,'Location','northeast','Interpreter','latex','FontSize',10,'Box','off');
set(ax3,'YGrid','on','Box','on','YLim',params.esf_ylim);
set(ax3,'FontName','Times New Roman','FontSize',10);

%% ----------------------
% 6. y 轴标签微调
% ----------------------
% exportgraphics 的紧裁切会在最左元素上固定裁掉约 4 pt，正好切掉 y 轴标签的字头。
% 这里在标签左侧 2% 轴宽处放一个白色（不可见）竖排字符，把裁切边界推到标签之外。
% 偏移量按标签实际位置计算，避免多留空白把图幅压小。
for ax = [ax1, ax3]
    yl = ax.YLabel;
    yl.Units = 'normalized';
    ylab_x = yl.Position(1);
    text(ax, ylab_x - 0.026, 0.5, '|', 'Units', 'normalized', ...
        'Color', 'w', 'FontSize', 11, 'Rotation', 90, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Clipping', 'off', 'Interpreter', 'none');
end

%% ----------------------
% 5. 完成（导出由调用方负责，避免在此处写死文件名）
% ----------------------
end