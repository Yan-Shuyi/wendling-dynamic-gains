function [spk_rate, spk_t, dbg] = detect_spikes(x, Fs, p, ref)
% DETECT_SPIKES  斜率 + 幅度 + 波峰因数 + 不应期 的棘波计数（与采样率无关）。
%
%   [spk_rate, spk_t, dbg] = detect_spikes(x, Fs, p, ref)
%
% 输入
%   x    一段 EEG（mV），内部自己去直流
%   Fs   采样率（Hz）；内部低通到 spk_lp 并统一重采样到 512 Hz
%   p    classify_params_v2()
%   ref  make_reference() 的输出（用 ref.slope_p99 / ref.rms）
%
% 四条判据（缺一条都会被实测反例打穿，依据都在 tests/results/）
%   1. 斜率 ≥ ref.slope_p99 × spk_slope_k          —— 与参考段同一把尺子
%      不低通直接算斜率时，256 Hz 白噪声的 p99 就有 766 mV/s，比 8 mV 棘波的最大
%      斜率（340）还大，检不出任何棘波（calib_spike.txt）。
%   2. 峰值 ≥ ref.rms × spk_amp_k                  —— 幅度门槛，压噪声尖刺
%   3. 斜率 ≥ 本窗斜率 p90 × spk_win_k             —— 压"纯正弦被数成棘波"；
%      没有它，4 Hz 大振幅节律会刷出 6–13 Hz 的假棘波（timeline.txt 10–45 s）
%   4. 峰值 ≥ 本窗 rms × spk_crest                 —— 波峰因数；
%      没有它，低幅噪声窗会刷出 0.2–1 Hz 的假棘波（timeline.txt 45–85 s）
[d, xb, Fs_e] = slope_series(x, Fs, p.spk_lp);

thr_slope = ref.slope_p99 * p.spk_slope_k;
thr_amp   = ref.rms * p.spk_amp_k;
thr_win   = prctile(d, 90) * p.spk_win_k;
thr_slope = max(thr_slope, thr_win);
thr_crest = p.spk_crest * sqrt(mean(xb.^2));

cand  = find(d > thr_slope);
spk_t = [];
min_gap = max(1, round(p.spk_refr * Fs_e));
half    = round(0.02 * Fs_e);              % ±20 ms 内找峰
last = -inf;
for k = 1:numel(cand)
    if cand(k) - last < min_gap
        continue;                          % 同一棘波的上升/下降沿只算一次
    end
    lo = max(1, cand(k) - half);
    hi = min(numel(xb), cand(k) + half);
    [pk, im] = max(xb(lo:hi));
    if pk < thr_amp || pk < thr_crest
        continue;                          % 幅度不够 / 不是孤立尖峰
    end
    spk_t(end+1) = (lo + im - 1)/Fs_e;     %#ok<AGROW>
    last = cand(k);
end

spk_rate = numel(spk_t) / (numel(xb)/Fs_e);
dbg = struct('thr_slope', thr_slope, 'thr_amp', thr_amp, ...
             'thr_win', thr_win, 'thr_crest', thr_crest, 'n_cand', numel(cand));
end
