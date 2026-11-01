function [lab, sub, f, score] = classify_v2(x, Fs, p, ref)
% CLASSIFY_V2  按 Wendling 2005 的四个期判一个窗；判据全部来自文献频带。
%
%   [lab, sub, f, score] = classify_v2(x, Fs, p, ref)
%
% 输入
%   x    一段 EEG（mV）
%   Fs   采样率（Hz）
%   p    classify_params_v2()
%   ref  make_reference() 的输出；缺省时用本窗自身统计（不推荐）
%
% 输出
%   lab     'intICTAL' | 'preONSET' | 'ONSET' | 'ICTAL' | 'postICTAL' | 'uncertain'
%   sub     亚类：ICTAL 分 'delta'(<4) / 'theta'(4–8) / 'alpha'(8–13)；
%           ONSET 分 'lvf'（纯快活动）/ 'fast+spike'；
%           preONSET 分 'sporadic' / 'repetitive'；
%           uncertain 分 'margin'（判据贴阈值）/ 'offband'（高幅但不在四期频带内）
%   f       extract_features_v2 的全部特征
%   score   各判据的取值，便于画边界敏感性
%
% 判定顺序（这条顺序本身就是方法的一部分，写进 Methods）
%   1. postICTAL（可选，默认关）：RMS 低于参考段且主频 < 4 Hz
%   2. ONSET：P(15–40)/P(1–90) 高且不低于 P(4–10) —— 原文的"快放电"
%      **不设幅度门槛**：原文的 ONSET 是低电压快活动（LVF），加幅度门槛会把 LVF 判掉
%   3. ICTAL：P(4–10) 占优、主频不高于 13 Hz，且窗口 RMS ≥ 参考段 × thr_ictal_rms
%      **必须设幅度门槛**：不加的话低幅背景噪声（宽带、4–10 Hz 占比偶然过 40%）
%      会被判成 ICTAL —— 实测 baseline 窗有 8/21 被误判（见 tests/results/timeline.txt）
%   4. preONSET：以上都不满足但检出棘波
%   5. uncertain：高幅但不落在任何一期频带内 —— 诚实留白，不硬贴标签
%   6. intICTAL：其余
%   任一判据落在阈值 ±p.uncertain_band 内 → 'uncertain'，不硬贴标签。
if nargin < 4 || isempty(ref)
    ref = make_reference(x, Fs, p);
end

f = extract_features_v2(x, Fs, p);
xf = x(:) - mean(x(:));
[spk_rate, spk_t] = detect_spikes(xf, Fs, p, ref);
f.spk_rate = spk_rate;
f.spk_t = spk_t;

rr = f.rms / max(ref.rms, eps);
lab = 'intICTAL'; sub = '';
score = struct('rel_fast', f.rel_fast, 'rel_slow', f.rel_slow, ...
               'spk_rate', spk_rate, 'f_peak', f.f_peak, 'rms_ratio', rr);
ub = p.uncertain_band;

%% 1. postICTAL（默认关闭）
if p.post_enable && rr < p.post_rms_ratio && f.f_peak < p.post_f_peak_max
    lab = 'postICTAL'; return
end

%% 2. ONSET：快活动
if f.rel_fast >= p.thr_fast*(1-ub) && f.rel_fast >= f.rel_slow
    if f.rel_fast < p.thr_fast*(1+ub)
        lab = 'uncertain'; sub = 'margin'; return
    end
    lab = 'ONSET';
    if spk_rate >= p.thr_spike_lo, sub = 'fast+spike'; else, sub = 'lvf'; end
    return
end

%% 3. preONSET：棘波（在 ICTAL 之前判）
%%   理由：离散棘波和高幅节律的频谱会重叠（棘波串的 4–10 Hz 占比也能过 0.30），
%%   但"数得出棘波"这件事本身是排他的。先判棘波，2 Hz 棘波串就不会被当成 ICTAL 节律。
if spk_rate >= p.thr_spike_lo
    lab = 'preONSET';
    if spk_rate >= p.thr_spike_hi, sub = 'repetitive'; else, sub = 'sporadic'; end
    return
end

%% 4. ICTAL：4–10 Hz 节律 + 幅度门槛
if rr >= p.thr_ictal_rms && f.rel_slow >= p.thr_slow*(1-ub) && f.f_peak <= p.slow_peak_max
    if f.rel_slow < p.thr_slow*(1+ub) || rr < p.thr_ictal_rms*(1+ub)
        lab = 'uncertain'; sub = 'margin'; return
    end
    lab = 'ICTAL';
    if f.f_peak < 4
        sub = 'delta';        % 原文只写到 4–10 Hz；delta 是"继续变慢"的延伸，Methods 里说明
    elseif f.f_peak < 8
        sub = 'theta';
    else
        sub = 'alpha';
    end
    return
end

%% 5. 高幅但不在四期频带内 → 留白，不硬贴
if rr >= p.thr_ictal_rms
    lab = 'uncertain'; sub = 'offband'; return
end

%% 6. intICTAL
lab = 'intICTAL';
end
