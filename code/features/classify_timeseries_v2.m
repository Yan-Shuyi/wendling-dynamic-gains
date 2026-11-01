function [t_center, labs, idx_full, feats] = classify_timeseries_v2(x, Fs, p, ref)
% CLASSIFY_TIMESERIES_V2  用 v2 判据把整条 EEG 逐窗分类（Fig.4 的中层色条）。
%
%   [t_center, labs, idx_full, feats] = classify_timeseries_v2(x, Fs, p, ref)
%
% 输入
%   x     整段 EEG（1 × N，mV）
%   Fs    采样率（Hz）
%   p     classify_params_v2()
%   ref   参考段（make_reference 的输出）；省略时用整条轨迹前 10 s
%
% 输出
%   t_center  每个窗的中心时刻（s）
%   labs      每个窗的标签（cell）
%   idx_full  逐采样点的类别编号（1..numel(p.names)），NaN 表示窗未覆盖
%   feats     每窗的特征结构体数组（rms / f_peak / rel_slow / rel_fast / spk_rate）
%
% 做什么
%   10 s 窗、2 s 步长（与全文一致），窗内直接调 classify_v2。窗不越界、不做镜像填充：
%   边界处不足一个窗的部分留 NaN，画图时按空白处理，避免"填充出来的假标签"。
%
% 备注
%   这是 classify_eeg_timeseries 的 v2 版本。旧版按 5 s 窗 + 六类振幅分区阈值工作，
%   与 §3.5 的频带判据不是一套，已停用（保留文件仅为追溯旧图）。
if nargin < 4 || isempty(ref)
    k = min(numel(x), round(10*Fs));
    ref = make_reference(x(1:k), Fs, p);
end
x = x(:)';
win = p.win_sec; step = 2;
N = numel(x);
t0v = 0:step:(N/Fs - win);
nw = numel(t0v);

t_center = t0v + win/2;
labs = cell(1, nw);
feats = repmat(struct('rms', [], 'f_peak', [], 'rel_slow', [], 'rel_fast', [], 'spk_rate', []), 1, nw);
for w = 1:nw
    k1 = round(t0v(w)*Fs) + 1;
    k2 = k1 + round(win*Fs) - 1;
    [lab, ~, f] = classify_v2(x(k1:k2), Fs, p, ref);
    labs{w} = lab;
    feats(w).rms = f.rms;  feats(w).f_peak = f.f_peak;
    feats(w).rel_slow = f.rel_slow;  feats(w).rel_fast = f.rel_fast;
    feats(w).spk_rate = f.spk_rate;
end

idx_full = NaN(1, N);
for w = 1:nw
    k1 = round(t0v(w)*Fs) + 1;
    k2 = min(N, k1 + round(win*Fs) - 1);
    li = find(strcmp(p.names, labs{w}), 1);
    if isempty(li), li = numel(p.names); end
    idx_full(k1:k2) = li;
end
end
