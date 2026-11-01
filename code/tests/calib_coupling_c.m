function calib_coupling_c()
% CALIB_COUPLING_C  用"驱动—响应"曲线标定三节点网络要用的全局耦合常数 C。
%
% 背景：耦合项写成 C·W·S(y_out)，量纲是脉冲密度（pulses/s），与噪声均值 90 同量纲。
% 因此 C 不能凭直觉取 1；必须让"一条满权重连接"的贡献落在某个可解释的位置上。
% probe_coupling_scale 已给出三个候选参数点的 S(y_out)：
%   NIZ  S_mean = 0.058   （y_out 长期为负，发放率几乎为零）
%   PZ   S_mean = 0.185
%   SOZ  S_mean = 1.110   （放电时 p95 = 4.87，接近饱和）
%
% 本脚本扫 C，量 receiver 被 driver 拉动的三件事：
%   rms 比、主频、PLV，外加 v2 判据给出的标签构成——最后一项才是"是否被募集"的判据，
%   因为相位锁定（PLV）不等于被募集（见输出）。
%
% 结构：两节点单向，1 = driver（SOZ，静态增益，持续 4 Hz 节律），2 = receiver。
% 静态增益是为了把"耦合拉动"与"receiver 自己的增益级联"分开；动态复核另做。
%
% 输出：tests/results/calib_c_sweep.txt/.mat + figures/calib_coupling_c.png

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');
figdir = fullfile(root, '..', 'figures');

Cv = [0 5 7.5 10 15 20 30 40 60 80 120];
recv = {'PZ', 4.00, 35, 10; 'NIZ', 3.25, 45, 25};
DUR = 60; T0 = 20;

p = classify_params_v2();
fid = fopen(fullfile(resdir, 'calib_c_sweep.txt'), 'w');
fprintf(fid, '# driver = SOZ (A=4.60, B=25, G=10), static gains, W(2,1)=1, delay=0\n');
fprintf(fid, '# recv  C     rms    rms_ratio  f_peak  R_slow  PLV_dr  label        frac\n');

R = struct();
for k = 1:size(recv, 1)
    name = recv{k,1};
    refk = [];
    for ic = 1:numel(Cv)
        C = Cv(ic);
        cfg = default_config();
        cfg.node_count = 2;
        cfg.A = [4.60, recv{k,2}];
        cfg.B = [25,   recv{k,3}];
        cfg.G = [10,   recv{k,4}];
        cfg.model_dynamic = 0;         % 静态增益：只看耦合拉动
        cfg.esf_enable = false;
        cfg.simulation_time = DUR;
        cfg.resample_dy = false;
        cfg.coupling_gain = C;
        cfg.coupling_clip_negative = false;   % 负向照样进 S（S 单调，值为小）
        cfg.sc_matrix = [0 0; 1 0];           % 1 → 2，单向
        cfg.seed = 20261031;
        y = simulate_network(cfg);     % 2 × n
        fs = cfg.eeg_fs;

        yd = y(1, :); yr = y(2, :);
        if C == 0
            refk = make_reference(yr(1:round(10*fs)), fs, p);   % 该 receiver 自己的静息基线
        end

        t  = (0:numel(yr)-1)/fs;
        keep = t > T0;
        [fr, pf, rf, rs] = describe(yr(keep), fs);
        plv = plv_pair(yd(keep), yr(keep));
        [tc, labs] = classify_timeseries_v2(yr, fs, p, refk);
        sel = tc > T0;
        labs = labs(sel);
        [ul, ~, iu] = unique(labs);
        cnt = accumarray(iu(:), 1);
        [~, im] = max(cnt);
        majlab = ul{im}; majfrac = cnt(im)/numel(labs);

        R.(name).rms(ic) = sqrt(mean(yr(keep).^2));
        R.(name).fr(ic)  = fr;
        R.(name).rs(ic)  = rs;
        R.(name).plv(ic) = plv;
        R.(name).lab{ic} = majlab;
        R.(name).frac(ic) = majfrac;

        fprintf(fid, '%-5s %5.1f %7.3f %9.3f %7.2f %7.3f %7.3f  %-11s %5.2f\n', ...
            name, C, R.(name).rms(ic), R.(name).rms(ic)/R.(name).rms(1), fr, rs, plv, majlab, majfrac);
    end
end
fclose(fid);
R.Cv = Cv; R.recv = recv; R.dur = DUR;
save(fullfile(resdir, 'calib_c_sweep.mat'), 'R');

%% ---------------- 判据：PZ 被募集、NIZ 不被募集 ----------------
% 募集 = 多数窗判为 preONSET 或 ICTAL，且占多数窗比例 ≥ 0.6
isRecruited = @(lab, frac) any(strcmp(lab, {'preONSET', 'ICTAL'})) && frac >= 0.6;
okPZ  = false(1, numel(Cv)); okNIZ = false(1, numel(Cv));
for ic = 1:numel(Cv)
    okPZ(ic)  = isRecruited(R.PZ.lab{ic},  R.PZ.frac(ic));
    okNIZ(ic) = isRecruited(R.NIZ.lab{ic}, R.NIZ.frac(ic));
end
idxPZ = find(okPZ); idxNIZ = find(okNIZ);

fprintf('\n--- 扫 C 结果（driver = SOZ 静态，W=1）---\n');
fprintf('%-6s %8s %8s %8s %8s %8s  %-10s %-10s\n', 'C', 'rmsPZ', 'plvPZ', 'rmsNIZ', 'plvNIZ', 'rslowPZ', 'labPZ', 'labNIZ');
for ic = 1:numel(Cv)
    fprintf('%-6.1f %8.3f %8.3f %8.3f %8.3f %8.3f  %-10s %-10s\n', Cv(ic), ...
        R.PZ.rms(ic), R.PZ.plv(ic), R.NIZ.rms(ic), R.NIZ.plv(ic), R.PZ.rs(ic), ...
        R.PZ.lab{ic}, R.NIZ.lab{ic});
end
if isempty(idxPZ)
    fprintf('\nPZ: 在所扫 C 范围内未被募集\n');
else
    fprintf('\nPZ 被募集: C ∈ [%g, %g]\n', Cv(idxPZ(1)), Cv(idxPZ(end)));
end
if isempty(idxNIZ)
    fprintf('NIZ: 全程未被募集\n');
else
    fprintf('NIZ 被募集: C ∈ [%g, %g] —— 该区间内选择性必须由 W 承担\n', Cv(idxNIZ(1)), Cv(idxNIZ(end)));
end

%% ---------------- 等效输入尺度 ----------------
% 放电时 S_mean(SOZ) = 1.11；C·1·1.11 / 90 就是"一条满权重连接贡献了几个噪声均值"
sbar = 1.110;
fprintf('\n等效尺度：C·W·S_mean / noise_mean（W=1）\n');
for ic = 1:numel(Cv)
    fprintf('  C=%5.1f → %5.3f 个噪声均值\n', Cv(ic), Cv(ic)*sbar/90);
end

%% ---------------- 图 ----------------
plot_calib_c(R);

fprintf('\n已写 results/calib_c_sweep.txt/.mat 与 figures/calib_coupling_c.png\n');
end

% ======================================================================
function [fpeak, pf, rf, rs] = describe(y, fs)
% 描述一段波形：主频、谱峰占比、快/慢占比
y = y(:)' - mean(y);
n = numel(y);
NFFT = 2^nextpow2(n);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));      % Hann 窗，手写以免依赖工具箱
Y = fft(y(:) .* w, NFFT);
P = abs(Y(1:NFFT/2+1)).^2;
fr = (0:NFFT/2)' * fs / NFFT;
band = @(a, b) trapz(fr(fr>=a & fr<=b), P(fr>=a & fr<=b));
Ptot = band(1, 90);
[pk, ip] = max(P(fr >= 1 & fr <= 60));
frr = fr(fr >= 1 & fr <= 60);
fpeak = frr(ip);
pf = pk / Ptot;
rf = band(15, 40) / Ptot;
rs = band(4, 10) / Ptot;
end

function v = plv_pair(a, b)
% Hilbert 相位差的 PLV
a = a(:)' - mean(a); b = b(:)' - mean(b);
ha = hilbert(a); hb = hilbert(b);
d = angle(ha) - angle(hb);
v = abs(mean(exp(1i * d)));
end
