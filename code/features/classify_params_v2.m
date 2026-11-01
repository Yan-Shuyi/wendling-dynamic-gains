function p = classify_params_v2()
% CLASSIFY_PARAMS_V2  Wendling 2005 四期判据的可执行版本：阈值全部无量纲。
%
%   与 v1（classify_params.m）的区别
%     1. 判据来自文献频带，不来自手工试出来的 mV 数：
%          ONSET  = 快放电，15–40 Hz（原文 beta / low-gamma）
%          ICTAL  = 变慢后的节律，4–10 Hz（原文 theta / alpha）
%          preONSET = 高幅棘波出现
%          intICTAL = 以上都无
%     2. 幅度类判据一律先除以参考段（同一仿真的发作间期基线）→ 换噪声、
%        换增益律、换 k 都不用重调阈值（v1 的 C8）。
%     3. 棘波用**斜率**（mV/s）判定，不用逐点差分 → 与采样率无关（v1 的 C5）。
%     4. 先减窗内均值再算所有特征；直流偏置单独输出为 depolarization，
%        不参与分期（v1 的 C1/C2）。
%
%   窗口
%     Wendling 2005 用的是 10 s 片段，这里默认也 10 s（v1 是 5 s）。
%
%   用法
%     p = classify_params_v2();
%     ref = make_reference(x_baseline, Fs, p);      % 参考段统计量
%     [lab, sub, f] = classify_v2(x_win, Fs, p, ref);
%
p.win_sec   = 10;      % 分析窗长（s），对齐 Wendling 2005 的 10 s 片段
p.step_sec  = 1;       % 滑窗步长（s）
p.hp_cut    = 0.5;     % 去直流用的高通截止（Hz）；0 表示只减均值
                        % 不能取 1 Hz：1 Hz 本身就在 delta 边上，会把 1 Hz 节律整段滤掉

%% 频带定义（Hz）——直接取自 Wendling 2005 原文
p.bands = { ...
    'slow',   1,  4; ...    % delta
    'theta',  4,  8; ...
    'alpha',  8, 13; ...
    'beta',  13, 30; ...
    'gamma', 30, 40};       % 原文 low-gamma 上界取 40
p.fast_band = [15 40];      % ONSET 判据：原文 "15 to 40 Hz"
p.slow_band = [4 10];       % ICTAL 判据：原文 "4 to 10 Hz"

%% 判据阈值（全部无量纲）
p.thr_fast   = 0.35;   % P(15–40)/P(1–90) 超过此值 → 快活动（标定：发作间期应 <0.15）
p.thr_slow   = 0.30;   % P(4–10)/P(1–90) 超过此值且主频不高于 13 Hz → 节律期
p.thr_dom    = 1.5;    % 主频带功率 / 次高带功率 超过此值才算"主导"（避免边界抖动）
p.thr_ictal_rms = 2.5; % ICTAL 的幅度门槛：窗 RMS ≥ 参考段 RMS × 此值
                       % 若不设，低幅宽带背景（4–10 Hz 占比偶然 >40%）会被判成 ICTAL
p.slow_peak_max = 13;  % ICTAL 的 f_peak 上界（Hz），再高就算快活动
p.thr_ictal_crest = 2.5; % ICTAL 还要求波峰因数 ≤ 此值：节律是"连续的鼓"（≈1.4），
                        % 棘波串是"孤立的尖"（≥3）——不加这条，2 Hz 棘波串会被判成 ICTAL

%% 棘波判据（与采样率无关）
p.spk_lp        = 40;   % 算斜率前先低通到此频率（棘波是 40 Hz 以下的形态）→ 与 Fs 无关
p.spk_slope_k   = 1.8;  % 斜率阈值 = ref.slope_p99 × k（mV/s）；校准见 tests/results/calib_spike.txt
                        % 取 1.8 的理由：8 mV 棘波是白噪声参考段 p99 的 2.3 倍、1/f 参考段的 20 倍，
                        % 取 3 时白噪声参考段检不出；1.8 时两种参考段都检得出且不误检
p.spk_amp_k     = 3;    % 且峰值 ≥ ref.rms × 此值，防噪声尖刺误检
p.spk_win_k     = 1.6;  % 且斜率 ≥ 本窗自身斜率 p90 × 此值 → 排除"纯正弦被数成棘波"
p.spk_crest     = 3;    % 且峰值 ≥ 本窗 rms × 此值（波峰因数）→ 排除低幅噪声刷假棘波
p.spk_refr      = 0.05; % 不应期（s），50 ms 内不重复计数
p.thr_spike_lo  = 0.2;  % 棘波率 ≥ 此值（Hz）→ preONSET（0.2 Hz ≈ 每 5 s 一个）
p.thr_spike_hi  = 0.5;  % 棘波率 ≥ 此值 → preONSET-repetitive（否则 sporadic）

%% postICTAL 判据（默认关闭：先要确认模型里真有"发作后抑制"）
p.post_enable = false;
p.post_rms_ratio = 0.7;  % 窗 RMS < 基线 RMS × 此值 → 抑制
p.post_f_peak_max = 4;   % 且主频 < 4 Hz

%% 不确定带：判据落在阈值 ±此比例内 → 'uncertain'（不要硬贴标签）
p.uncertain_band = 0.10;

%% 类别名（对外显示用；与 Wendling 2005 原文一致）
p.names = {'intICTAL','preONSET','ONSET','ICTAL','postICTAL','uncertain'};
end
