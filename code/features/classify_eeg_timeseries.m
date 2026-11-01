% ==============================
% plaintext: 长时程EEG逐窗分类主函数
% 支持滑窗步长，镜像填充边界，可选插值，输出类型序列
% ==============================
function varargout = classify_eeg_timeseries(output, Fs, params, do_plot)
% CLASSIFY_EEG_TIMESERIES  Sliding-window classification of a whole EEG trace.
%
%   [t_center, type_seq, type_full] = classify_eeg_timeseries(output, Fs, params, do_plot)
%
% 输入
%   output   整段 EEG（1 × N，mV）
%   Fs       采样率（Hz）
%   params   classify_params() 的阈值与滑窗设置
%   do_plot  是否顺手画一张带类型标注的图（默认 false；论文的图用
%            plot_pub_eeg_esf 单独画，所以这里是 false）
%
% 输出
%   t_center   每个滑窗的中心时刻（s），1 × K
%   type_seq   每个滑窗的类别编号，1 × K（1..6 或 NaN）
%   type_full  逐采样点的类别序列，1 × N；同一窗内的点取同一类别，
%              直接拿去给 plot_pub_eeg_esf 画中间的彩色类型条
%
% 做什么
%   1. 按 params.win_len / win_step 切滑窗，边界处做镜像填充；
%   2. 每个窗调用 classify_single 得到类别；
%   3. 把窗级类别展开成逐点序列 type_full。
%
% 调用关系
%   上游：run_fig2.m（生成 Fig. 4 中间那条类型条）
%   下游：classify_single、plot_timeseries_result（仅 do_plot = true 时）
%
% 备注
%   这是"增益演化"与"电生理表型"之间唯一的连接点：Fig. 2 的树在这里被真正
%   用上，Fig. 4 中间的类型条就是它的输出。
if nargin < 4
    do_plot = false;
end
output = output(:);
N = length(output);

%% ----------------------
% 1. 滑窗参数计算
% ----------------------
win_half_samples = round(params.win_half * Fs);
step_samples = round(params.step_sec * Fs);
step_samples = max(step_samples, 1); % 兜底：步长不小于1采样点

% 镜像填充边界（保证首尾窗口完整）
pad_front = output(win_half_samples:-1:1);
pad_back = output(end:-1:end-win_half_samples+1);
output_pad = [pad_front; output; pad_back];

%% ----------------------
% 2. 生成滑窗中心索引
% ----------------------
all_center_idx = (win_half_samples+1):(win_half_samples+N);
all_t_center = (0:N-1)/Fs;
selected_center_idx = all_center_idx(1:step_samples:end);
t_center = all_t_center(1:step_samples:end);
n_wins = length(selected_center_idx);

if n_wins == 0
    error('步长过大，无法生成有效窗口！请减小params.step_sec');
end

%% ----------------------
% 3. 逐窗分类
% ----------------------
type_seq = zeros(n_wins, 1);
fprintf('开始逐窗分类：共%d个窗口（步长%.2fs，计算量约为逐点模式的%.4f%%）\n', ...
    n_wins, params.step_sec, 100*n_wins/N);

for k = 1:n_wins
    % 提取当前窗口
    idx_start = selected_center_idx(k) - win_half_samples;
    idx_end = selected_center_idx(k) + win_half_samples;
    seg = output_pad(idx_start:idx_end);
    
    % 提取特征+分类
    features = extract_features(seg, params);
    type_seq(k) = classify_single(features, params);
    
    % 进度提示
    if mod(k, 100) == 0
        fprintf('进度：%d/%d (%.1f%%)\n', k, n_wins, 100*k/n_wins);
    end
end
fprintf('逐窗分类完成。\n');

%% ----------------------
% 4. 可选：插值到原始信号长度
% ----------------------
type_full = [];
if params.interp_to_signal
    fprintf('将类型序列插值到原始信号长度...\n');
    type_full = interp1(t_center, type_seq, all_t_center, 'nearest', 'extrap');
    type_full = type_full(:);
end

%% ----------------------
% 5. 输出结果
% ----------------------
varargout{1} = t_center;
varargout{2} = type_seq;
if nargout >= 3
    varargout{3} = type_full;
end

%% ----------------------
% 6. 可选绘图
% ----------------------
if do_plot
    if ~isempty(type_full)
        plot_timeseries_result(output, Fs, all_t_center, type_full, params);
    else
        plot_timeseries_result(output, Fs, t_center, type_seq, params);
    end
end
end