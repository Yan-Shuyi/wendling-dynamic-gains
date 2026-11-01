% ==============================
% plaintext: 特征提取逻辑
% 1. 基础时域特征：min_val（信号最小值），用于一级分类、RMS、线尺度
% 2. 棘波特征：频谱峰值数量
% 3. 频域特征：
%    - 主频、高频(15-90Hz)能量占比（区分Type5）
%    - 频谱显著峰数量（区分Type2/3，过滤噪声峰）
%    - 频谱90%能量带宽（区分Type4/6，纯净正弦带宽极窄）
% 4. 所有特征均为无量纲或物理量，无归一化，方便阈值校准
% ==============================
% ==============================
% plaintext: 特征提取（严格对齐分类规则所需特征）
% 提取所有分类必需的特征，无冗余，无遗漏
% ==============================
function features = extract_features(signal, params)
% EXTRACT_FEATURES  Compute the time/frequency features used by the classifier.
%
%   feat = extract_features(x, Fs, params)
%
% 输入
%   x       一个滑窗内的 EEG 片段（1 × n，mV）
%   Fs      采样率（Hz）
%   params  classify_params() 的结构体（用到去趋势、频谱估计等设置）
%
% 输出
%   feat    结构体，逐项特征：
%             mean_val / min_val / max_val / range_val   时域幅值统计
%             rms                                         均方根
%             line_length   线尺度 L = mean(|diff(x)|)，即李涵论文图 3-24 用的 L
%             dom_freq      主频（Hz）
%             peak_num      显著谱峰个数（已过滤噪声小峰）
%             hf_ratio      15 Hz 以上能量占比
%             bw90          90% 能量带宽（Hz）
%
% 做什么
%   把一段 EEG 变成决策树需要的那几个数。全部特征都是物理量或无量纲量，
%   不做归一化，方便直接对照阈值。
%
% 调用关系
%   上游：classify_single（每个滑窗调一次）
%   下游：无
%
% 备注
%   line_length 已经算出来了，但 classify_single.m 目前没用它——李涵论文的
%   树是用 L 分左右枝的。这个"算了不用"的状态见 log 的 D2，等你决定是否接入。
signal = signal(:);
N = length(signal);

%% ----------------------
% 1. 时域特征（低幅组核心特征）
% ----------------------
features.mean_val = mean(signal);                  % 均值（一级分类核心）
features.min_val = min(signal);                    % 最小值（仅在低幅组使用，符合“负值仅存于123”规则）
features.max_val = max(signal);                    % 最大值
features.range_val = features.max_val - features.min_val; % 极差（区分Type1和2/3）
features.rms = rms(signal);                        % 均方根幅值
features.line_length = sum(abs(diff(signal)))/N;   % 线尺度（波形复杂度）

% 差分特征（检测棘波陡峭程度）
diff_sig = abs(diff(signal));
features.diff_max = max(diff_sig);                 % 差分最大值
features.diff_std = std(diff_sig);                 % 差分标准差（用于自适应阈值）

%% ----------------------
% 2. 频域特征（高幅组核心特征）
% ----------------------
window_size = min(1024, N);
overlap = floor(window_size * 0.5);
nfft = max(1024, 2^nextpow2(window_size));
[pxx, f] = pwelch(signal, hanning(window_size), overlap, nfft, params.Fs);
idx_band = f >= 1 & f <= 90;
f = f(idx_band);
pxx = pxx(idx_band);
pxx_norm = pxx / sum(pxx); % 归一化到总能量为1

% 主频
[~, peak_idx] = max(pxx_norm);
features.f_peak = f(peak_idx);

% 高频(15-90Hz)能量占比（区分Type5）
idx_high = f >= params.high_freq_start & f <= 90;
features.high_freq_ratio = sum(pxx_norm(idx_high));

% 频谱90%能量带宽（区分Type4/6）
cum_energy = cumsum(pxx_norm);
idx_90 = find(cum_energy >= 0.9, 1, 'first');
if ~isempty(idx_90)
    features.spec_bw = f(idx_90) - f(1);
else
    features.spec_bw = f(end) - f(1); % 找不到90%点时用全带宽兜底
end

% 频谱显著尖峰数量（区分Type2/3）
min_peak_height = params.min_peak_height_ratio * max(pxx_norm);
[~, peak_locs] = findpeaks(pxx_norm, f, ...
    'MinPeakHeight', min_peak_height, ...
    'MinPeakDistance', params.min_peak_dist); % 用频率轴，单位为Hz
features.peak_num = length(peak_locs);
end