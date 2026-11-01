function calib_coupling_seed()
% CALIB_COUPLING_SEED  C 的噪声稳健性：同一 C 下换 10 个随机种子，看标签会不会翻。
%
% 为什么必须做：calib_coupling_c 的每一格都是单个 60 s 实现。若某个 C 下
% "PZ 被募集"只在特定噪声实现里成立，那这个 C 就不能定下来。
%
% 结构同 calib_coupling_c（单向、静态增益、W=1），C 取 15/20/25/30/40，种子 10 个。
% 判定"被募集"用两个条件同时成立：标签 ∈ {preONSET, ICTAL} 且 rms 比 ≥ 1.5。
%
% 输出：tests/results/calib_c_seed.txt/.mat

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

Cv = [15 20 25 30 40];
seeds = 20261031 + (0:9);
recv = {'PZ', 4.00, 35, 10; 'NIZ', 3.25, 45, 25};
DUR = 60; T0 = 20;

p = classify_params_v2();
fid = fopen(fullfile(resdir, 'calib_c_seed.txt'), 'w');
fprintf(fid, '# C  seed  recv   rms   ratio  f_peak  R_slow  label      frac  recruited\n');

S = struct();
for k = 1:size(recv, 1)
    name = recv{k,1};
    for ic = 1:numel(Cv)
        for is = 1:numel(seeds)
            cfg = default_config();
            cfg.node_count = 2;
            cfg.A = [4.60, recv{k,2}];
            cfg.B = [25,   recv{k,3}];
            cfg.G = [10,   recv{k,4}];
            cfg.model_dynamic = 0;
            cfg.esf_enable = false;
            cfg.simulation_time = DUR;
            cfg.resample_dy = false;
            cfg.coupling_gain = Cv(ic);
            cfg.coupling_clip_negative = false;
            cfg.sc_matrix = [0 0; 1 0];
            cfg.seed = seeds(is);
            y = simulate_network(cfg);
            fs = cfg.eeg_fs;

            yr = y(2, :);
            % 参考段仍取该 receiver 在 C=0、同一 seed 下的前 10 s
            cfg0 = cfg; cfg0.coupling_gain = 0;
            y0 = simulate_network(cfg0);
            ref = make_reference(y0(2, 1:round(10*fs)), fs, p);

            t = (0:numel(yr)-1)/fs;  keep = t > T0;
            [fr, ~, ~, rs] = describe(yr(keep), fs);
            [tc, labs] = classify_timeseries_v2(yr, fs, p, ref);
            labs = labs(tc > T0);
            [ul, ~, iu] = unique(labs); cnt = accumarray(iu(:), 1);
            [~, im] = max(cnt); majlab = ul{im}; majfrac = cnt(im)/numel(labs);

            ratio = sqrt(mean(yr(keep).^2)) / sqrt(mean(y0(2, keep).^2));
            recruited = any(strcmp(majlab, {'preONSET','ICTAL'})) && ratio >= 1.5;

            S.(name).ratio(ic, is) = ratio;
            S.(name).fr(ic, is)    = fr;
            S.(name).rs(ic, is)    = rs;
            S.(name).rec(ic, is)   = recruited;
            S.(name).lab{ic, is}   = majlab;

            fprintf(fid, '%4g %10d %-5s %6.3f %6.2f %6.2f %6.3f  %-10s %4.2f  %d\n', ...
                Cv(ic), seeds(is), name, sqrt(mean(yr(keep).^2)), ratio, fr, rs, majlab, majfrac, recruited);
        end
    end
end
fclose(fid);
S.Cv = Cv; S.seeds = seeds; S.recv = recv;
save(fullfile(resdir, 'calib_c_seed.mat'), 'S');

fprintf('\n--- 每 C 的募集率（10 个种子）---\n');
fprintf('%-5s %10s %10s   %-22s %-22s\n', 'C', 'PZ_rate', 'NIZ_rate', 'PZ labels', 'NIZ labels');
for ic = 1:numel(Cv)
    lp = S.PZ.lab(ic, :); ln = S.NIZ.lab(ic, :);
    fprintf('%-5g %9.0f%% %9.0f%%   %-22s %-22s\n', Cv(ic), ...
        100*mean(S.PZ.rec(ic,:)), 100*mean(S.NIZ.rec(ic,:)), ...
        strjoin(unique(lp), ','), strjoin(unique(ln), ','));
end
fprintf('\n--- rms 比（中位 [最小, 最大]）---\n');
for ic = 1:numel(Cv)
    fprintf('C=%-4g PZ %5.2f [%5.2f,%5.2f]   NIZ %5.2f [%5.2f,%5.2f]\n', Cv(ic), ...
        median(S.PZ.ratio(ic,:)), min(S.PZ.ratio(ic,:)), max(S.PZ.ratio(ic,:)), ...
        median(S.NIZ.ratio(ic,:)), min(S.NIZ.ratio(ic,:)), max(S.NIZ.ratio(ic,:)));
end
fprintf('\n已写 results/calib_c_seed.txt/.mat\n');
end

function [fpeak, pf, rf, rs] = describe(y, fs)
y = y(:)' - mean(y);
n = numel(y);
NFFT = 2^nextpow2(n);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
Y = fft(y(:) .* w, NFFT);
P = abs(Y(1:NFFT/2+1)).^2;
fr = (0:NFFT/2)' * fs / NFFT;
band = @(a, b) trapz(fr(fr>=a & fr<=b), P(fr>=a & fr<=b));
Ptot = band(1, 90);
[pk, ip] = max(P(fr >= 1 & fr <= 60));
frr = fr(fr >= 1 & fr <= 60);
fpeak = frr(ip); pf = pk/Ptot;
rf = band(15,40)/Ptot; rs = band(4,10)/Ptot;
end
