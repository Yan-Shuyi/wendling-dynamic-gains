function test_classifier_boundaries()
% TEST_CLASSIFIER_BOUNDARIES  分期判据的三类检验：判据边界、噪声鲁棒性、独立描述子。
%
%   运行：run_checks   （或 addpath('C:/Users/26442/WorkBuddy/Claw/code'); setup_paths; test_classifier_boundaries）
%
% 为什么需要这三项（论文 Methods 要写"分类器已验证"，就必须有可复现的证据）
%   1. 原型检验（test_classifier_archetypes）只说明"极端情况判得对"，不说明判据在
%      **阈值附近**的行为。这里让某个判据量连续穿过阈值，看标签是否只在声称的
%      不确定带内变化，且过渡单调、不来回跳。
%   2. 判据的输入是带噪信号。给同一批原型加不同强度的噪声，统计标签保持率，
%      论文里才能写"标签在某 SNR 以上稳定"，这比"阈值取自文献"有力。
%   3. 最硬的一条：用一个**独立于判据的特征族**（自相关主周期）去核对同一批仿真，
%      两套毫不相干的量给出一致判断，才说明标签不是频带比值的凑数。
%
% 另外顺手量了一条**适用域**：多窄的瞬变会被误当成快活动（见 [A4]）。
%   这一条不写入"正确性"，写成 Methods 里的限制说明。
%
% 输出：results/classifier_boundaries.txt
here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
setup_paths();
fid = fopen(fullfile(outdir, 'classifier_boundaries.txt'), 'w');
p  = classify_params_v2();
Fs = 256;
T  = p.win_sec;
t  = (0:1/Fs:T-1/Fs)';

%% 参考段：低幅背景（与原型检验同一构造，保证口径一致）
rng(7);
base = 0.8*randn(round(T*Fs),1) + 0.5*sin(2*pi*1*t);
ref  = make_reference(base, Fs, p);
ub   = p.uncertain_band;

fprintf(fid, '分期判据边界 / 鲁棒性 / 独立描述子检验（win = %.0f s，Fs = %d）\n', T, Fs);
fprintf(fid, '参考段：rms %.3f mV，slope_p99 %.1f mV/s；不确定带 ±%.0f%%\n', ...
    ref.rms, ref.slope_p99, 100*ub);

%% ============================================================
%% [A1] 快活动判据 thr_fast：让 P(15-40)/P(1-90) 连续穿过 0.35
%%      构造 = 固定宽带背景 + 逐步增强的 25 Hz 快波（背景不当 6 Hz 正弦，
%%      否则 relSlow 一直压着 relFast，ONSET 分支永远进不去）
%% ============================================================
fprintf(fid, '\n[A1] 快活动判据 thr_fast = %.2f（不确定带 %.3f–%.3f）\n', ...
    p.thr_fast, p.thr_fast*(1-ub), p.thr_fast*(1+ub));
fprintf(fid, '%8s %9s %9s %9s %8s  %-16s\n', 'a25', 'relFast', 'relSlow', 'fPeak', 'crest', 'label');
rng(21);
bg = 0.8*randn(size(t));
flip_fast = []; prev = '';
for a = [0 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.7 0.8 0.9 1 1.2 2 4 8]
    x = bg + a*sin(2*pi*25*t);
    [lab, sub, f] = classify_v2(x, Fs, p, ref);
    s = sub_label(lab, sub);
    fprintf(fid, '%8.1f %9.3f %9.3f %9.1f %8.2f  %-16s\n', a, f.rel_fast, f.rel_slow, f.f_peak, f.crest, s);
    if ~isempty(prev) && ~strcmp(s, prev), flip_fast(end+1) = f.rel_fast; end
    prev = s;
end
if isempty(flip_fast)
    fprintf(fid, '  过渡处的 relFast：无过渡（检查构造）\n');
else
    fprintf(fid, '  过渡处的 relFast：%s（应落在 %.3f–%.3f 内）\n', ...
        mat2str(round(flip_fast,3)), p.thr_fast*(1-ub), p.thr_fast*(1+ub));
end

%% ============================================================
%% [A2] 慢节律判据 thr_slow：让 P(4-10)/P(1-90) 穿过 0.30
%%      构造 = 可变幅度 6 Hz 慢波 + 固定背景噪声；噪声幅度调到让 rms 比全程 >2.5，
%%      这样过渡只由 relSlow 决定，而不是被幅度门槛抢先
%% ============================================================
fprintf(fid, '\n[A2] 慢节律判据 thr_slow = %.2f（不确定带 %.3f–%.3f），幅度门槛 rms >= %.1f × 参考\n', ...
    p.thr_slow, p.thr_slow*(1-ub), p.thr_slow*(1+ub), p.thr_ictal_rms);
fprintf(fid, '%8s %9s %9s %9s %8s  %-16s\n', 'a6', 'relSlow', 'relFast', 'rms比', 'fPeak', 'label');
% a6 = 6 Hz 慢波幅度（mV）；噪声标准差随 k 同步增大，保持总 rms 不变
rng(11);
flip_slow = []; prev = '';
% 构造要点：让总 rms 基本不变（≈4.2–4.6 mV，全程高于幅度门槛），
% 只让"6 Hz 慢波 / 宽带噪声"的比例变化，从而 relSlow 连续穿过 0.30。
for k = 0:0.05:1
    a6 = 6.5*(1-k) + 0.5;   sg = 1.0 + 3.2*k;
    x = 12 + a6*sin(2*pi*6*t) + sg*randn(size(t));
    [lab, sub, f] = classify_v2(x, Fs, p, ref);
    s = sub_label(lab, sub);
    fprintf(fid, '%8.1f %9.3f %9.3f %9.2f %8.1f  %-16s\n', ...
        a6, f.rel_slow, f.rel_fast, f.rms/ref.rms, f.f_peak, s);
    if ~isempty(prev) && ~strcmp(s, prev), flip_slow(end+1) = f.rel_slow; end
    prev = s;
end
if isempty(flip_slow)
    fprintf(fid, '  过渡处的 relSlow：无过渡（检查构造）\n');
else
    fprintf(fid, '  过渡处的 relSlow：%s（应落在 %.3f–%.3f 内）\n', ...
        mat2str(round(flip_slow,3)), p.thr_slow*(1-ub), p.thr_slow*(1+ub));
end

%% ============================================================
%% [A3] 棘波率判据：扫描真实棘波率（棘波宽度取生理值 20 ms）
%% ============================================================
fprintf(fid, '\n[A3] 棘波率判据 lo = %.2f Hz（preONSET）、hi = %.2f Hz（repetitive）\n', ...
    p.thr_spike_lo, p.thr_spike_hi);
fprintf(fid, '%8s %9s %9s %9s  %-18s\n', '真值Hz', '测得Hz', '误差', 'relFast', 'label');
for r = [0.05 0.10 0.15 0.19 0.21 0.30 0.45 0.55 0.70]
    x = spike_train(t, r, 0.02, 8);
    [lab, sub, f] = classify_v2(x, Fs, p, ref);
    fprintf(fid, '%8.2f %9.2f %9.2f %9.3f  %-18s\n', r, f.spk_rate, f.spk_rate - r, f.rel_fast, ...
        sub_label(lab, sub));
end

%% ============================================================
%% [A4] 适用域：棘波宽度多窄会被误判成快活动（写成限制，不写成正确性）
%% ============================================================
fprintf(fid, '\n[A4] 适用域：棘波宽度扫描（率固定 2 Hz）—— 过窄的瞬变频谱是宽带的，\n');
fprintf(fid, '     会被 ONSET 分支当成快活动；20 ms 以上的生理棘波不会。\n');
fprintf(fid, '%10s %10s %9s %8s  %-18s\n', '宽度ms', '半高宽ms', 'relFast', 'crest', 'label');
for w = [0.004 0.0067 0.01 0.02 0.04]
    x = spike_train(t, 2, w, 8);
    [lab, sub, f] = classify_v2(x, Fs, p, ref);
    fprintf(fid, '%10.1f %10.1f %9.3f %8.2f  %-18s\n', ...
        1000*w, 1000*2*sqrt(log(2))*w/2*2, f.rel_fast, f.crest, sub_label(lab, sub));
end

%% ============================================================
%% [B] 噪声鲁棒性：8 个原型 × 6 个噪声水平 × 5 次实现
%% ============================================================
fprintf(fid, '\n[B] 噪声鲁棒性（噪声标准差按信号 rms 的百分比加入；每档 5 次实现，表内为标签保持率）\n');
A = archetypes(t);
lvl = [0 25 50 100 200 400];
fprintf(fid, '%-28s', '原型');
for L = lvl, fprintf(fid, ' %5d%%', L); end
fprintf(fid, '\n%s\n', repmat('-', 1, 66));
keep_sum = zeros(1, numel(lvl));
for k = 1:size(A,1)
    nm = A{k,1}; x0 = A{k,2};
    fprintf(fid, '%-28s', nm);
    s0 = x0 - mean(x0); r0 = sqrt(mean(s0.^2));
    [lab0, sub0] = classify_v2(x0, Fs, p, ref);
    for L = 1:numel(lvl)
        nk = 0;
        for rep = 1:5
            rng(100*k + rep);
            x = x0 + lvl(L)/100*r0*randn(size(x0));
            [lab, sub] = classify_v2(x, Fs, p, ref);
            nk = nk + (strcmp(lab, lab0) && strcmp(sub, sub0));
        end
        fprintf(fid, ' %5.0f%%', 100*nk/5);
        keep_sum(L) = keep_sum(L) + nk/5;
    end
    fprintf(fid, '   (基准 %s)\n', sub_label(lab0, sub0));
end
fprintf(fid, '%-28s', '平均保持率');
for L = 1:numel(lvl), fprintf(fid, ' %5.0f%%', 100*keep_sum(L)/size(A,1)); end
fprintf(fid, '\n');

%% ============================================================
%% [C] 独立描述子交叉检验：自相关主周期 vs 判据标签
%%     参考段改用"低兴奋性参数的仿真输出"，与扫描脚本口径一致
%% ============================================================
fprintf(fid, '\n[C] 独立描述子交叉检验（自相关主周期：只用时域周期性，不用任何频带功率比值）\n');
refC = sim_reference(A0_default(), Fs);
fprintf(fid, '     参考段来自 A=4.60、B=50、G=25 的仿真：rms %.3f mV，slope_p99 %.1f mV/s\n', ...
    refC.rms, refC.slope_p99);
cells = { ...
     0  30    '快活动区'; ...
     5  45    '快活动区'; ...
    10  45    '快活动区'; ...
    25  10    '慢节律区'; ...
    30   5    '慢节律区'; ...
    25  25    '慢节律区'; ...
    45  20    '背景'; ...
    50  10    '背景'; ...
    35  10    '棘波区'};
fprintf(fid, '%4s %4s %-9s %8s %9s %9s %8s  %-18s %-9s %s\n', ...
    'B', 'G', '区域', 'rms', 'f_peak', 'f_acf', 'ACF峰', 'v2 标签', '独立判断', '一致');
nOk = 0; nTot = 0;
for k = 1:size(cells,1)
    cfg = default_config();
    cfg.node_count = 1; cfg.model_dynamic = 0; cfg.esf_enable = false;
    cfg.resample_dy = false; cfg.simulation_time = 45;
    cfg.A = A0_default(); cfg.B = cells{k,1}; cfg.G = cells{k,2};
    cfg.seed = 4242 + k;
    out = simulate_network(cfg);
    x = out(:)'; x = x(round(20*Fs)+1 : round(30*Fs));
    [lab, sub, f] = classify_v2(x, Fs, p, refC);
    [fA, aPeak] = acf_period(x, Fs);
    % 独立判据（只用自相关）：周期性显著（ACF 峰 >= 0.3）且主周期落在
    % δ–α（1–13 Hz）→ 慢节律；落在 β–γ（>13 Hz）→ 快活动；否则无节律。
    % 频段边界按生理划分，不沿用分类器自己的 4–10 / 15–40，避免自证。
    if aPeak >= 0.3 && fA >= 1 && fA <= 13
        indep = '慢节律';
    elseif aPeak >= 0.3 && fA > 13
        indep = '快活动';
    else
        indep = '无节律';
    end
    switch lab
        case 'ICTAL',    exp_indep = '慢节律';
        case 'ONSET',    exp_indep = '快活动';
        case 'intICTAL', exp_indep = '无节律';
        otherwise,       exp_indep = '';
    end
    if isempty(exp_indep)
        ok_s = 'n/a';
    else
        ok = strcmp(indep, exp_indep);
        nOk = nOk + ok; nTot = nTot + 1;
        ok_s = ternary(ok, '是', '否');
    end
    fprintf(fid, '%4.0f %4.0f %-9s %8.3f %9.1f %9.2f %8.2f  %-18s %-9s %s\n', ...
        cells{k,1}, cells{k,2}, cells{k,3}, f.rms, f.f_peak, fA, aPeak, ...
        sub_label(lab,sub), indep, ok_s);
end
fprintf(fid, '  一致 %d / %d（preONSET/uncertain 不参与一致性统计）\n', nOk, nTot);

fclose(fid);
fprintf('已写出：%s\n', fullfile(outdir, 'classifier_boundaries.txt'));
end

%% ================= 本地工具 =================
function a = A0_default()
a = 4.60;
end

function s = sub_label(lab, sub)
if isempty(sub), s = lab; else, s = [lab '-' sub]; end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function x = spike_train(t, rate, sigma, amp)
% SPIKE_TRAIN  规则间隔的高斯棘波串；sigma 是高斯宽度（s），生理棘波取 0.02
x = zeros(size(t));
for tc = 0.5 : 1/rate : t(end)
    x = x + amp*exp(-((t-tc)/sigma).^2);
end
end

function A = archetypes(t)
% ARCHETYPES  8 个"应该判成什么"已知的合成原型（与原型检验同一批，棘波宽 20 ms）
A = {};
A(end+1,:) = {'4 Hz 正弦 amp=8, dc=12',   12+8*sin(2*pi*4*t)};
A(end+1,:) = {'9 Hz 正弦 amp=5, dc=10',   10+5*sin(2*pi*9*t)};
A(end+1,:) = {'25 Hz 快波 amp=2, dc=0',   2*sin(2*pi*25*t)};
A(end+1,:) = {'35 Hz 快波 amp=2, dc=8',   8+2*sin(2*pi*35*t)};
A(end+1,:) = {'2 Hz 棘波串 amp=8',        spike_train(t,2,0.02,8)};
A(end+1,:) = {'0.3 Hz 棘波串 amp=8',      spike_train(t,0.3,0.02,8)};
rng(5);  A(end+1,:) = {'白噪声 sigma=0.8', 0.8*randn(size(t))};
A(end+1,:) = {'常数 10 mV',                10*ones(size(t))};
end

function ref = sim_reference(A0, Fs)
% SIM_REFERENCE  用低兴奋性、强抑制参数的仿真输出当参考段（"发作间期基线"）
cfg = default_config();
cfg.node_count = 1; cfg.model_dynamic = 0; cfg.esf_enable = false;
cfg.resample_dy = false; cfg.simulation_time = 40;
cfg.A = A0; cfg.B = 50; cfg.G = 25; cfg.seed = 99;
out = simulate_network(cfg);
x = out(:)';
x = x(round(20*Fs)+1 : round(30*Fs));
ref = make_reference(x, Fs, classify_params_v2());
end

function [f0, peak] = acf_period(x, Fs)
% ACF_PERIOD  自相关主周期：先找第一个负穿越，再在其后取第一个峰
%   f0 = 1/主周期（Hz）；peak = 该延迟处的自相关系数（节律性强度，0–1）
x = x(:) - mean(x(:));
[ac, lag] = xcorr(x, x, 'coeff');
keep = lag >= 0;
ac = ac(keep); lag = lag(keep);
k0 = find(ac < 0, 1, 'first');
if isempty(k0)                      % 没有负穿越 = 没有振荡，只有慢漂移/噪声
    f0 = NaN; peak = 0; return
end
hi = min(numel(ac), round(Fs/2));   % 最慢只看到 2 Hz
if k0 >= hi
    f0 = NaN; peak = 0; return
end
[pk, ip] = max(ac(k0:hi));
f0 = Fs / lag(k0 + ip - 1);
peak = pk;
end
