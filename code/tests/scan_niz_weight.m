function scan_niz_weight()
% SCAN_NIZ_WEIGHT  改变 SOZ→NIZ 的结构权重，定位"第三区保持安静"的上界。
%
%   scan_niz_weight()
%
% 目的
%   run_three_zone 用 W(3,1) = 0.15 时 NIZ 完全没被募集（Ge 一直是 3.25，
%   1 s rms 峰值 1.4 mV）。这个结果是"对的"，但只有当它落在一个稳定的窗口里
%   才说明问题；如果 0.2 就翻，那 0.15 只是侥幸。这里扫 0.15/0.30/0.50/0.80。
%
% 输出
%   tests/results/niz_weight.txt

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); setup_paths();
resdir = fullfile(root, 'tests', 'results');

T = 600; C = 15;
SOZ = [3.25 22  8 0.005];
PZ  = [3.25 32 12 0];
NIZ = [3.25 50 22 0];
D   = [0 0 0; 0.015 0 0; 0.025 0 0];

wlist = [0.15 0.30 0.50 0.80];

fid = fopen(fullfile(resdir, 'niz_weight.txt'), 'w');
fprintf(fid, '# C=%d, T=%d s, 只改 W(NIZ<-SOZ)\n', C, T);
fprintf(fid, '# w_niz   zone  rms_base  rms_max  Ge_max  t_Ge3.6  t_elec  ONSET_frac  ICTAL_frac\n');
fprintf('%-7s %-5s %9s %8s %8s %9s %8s %11s %10s\n', ...
    'w_niz','zone','rms_b','rms_max','Ge_max','t_Ge>3.6','t_elec','ONSET','ICTAL');

p = classify_params_v2();
nm = {'SOZ','PZ','NIZ'};
for wi = 1:numel(wlist)
    W = [0 0 0; 1 0 0; wlist(wi) 0 0];
    cfg = default_config();
    cfg.node_count = 3;
    cfg.A = [SOZ(1) PZ(1) NIZ(1)]; cfg.B = [SOZ(2) PZ(2) NIZ(2)]; cfg.G = [SOZ(3) PZ(3) NIZ(3)];
    cfg.gain_bias = [SOZ(4) PZ(4) NIZ(4)];
    cfg.model_dynamic = 1; cfg.gain_law = 'legacy'; cfg.gain_reset_tol = 0.05;
    cfg.esf_enable = false; cfg.simulation_time = T; cfg.resample_dy = true;
    cfg.seed = 20261101; cfg.coupling_gain = C;
    cfg.sc_matrix = W; cfg.delay_matrix = D;

    [y, dy, cfg] = simulate_network(cfg);
    fs = cfg.eeg_fs;
    w = round(1*fs); nr = floor(size(y,2)/w);
    rr = zeros(3, nr);
    for i = 1:nr
        seg = y(:,(i-1)*w+(1:w));
        rr(:,i) = sqrt(mean(seg.^2, 2));
    end
    tr = (0:nr-1) + 0.5;

    for i = 1:3
        Ge = dy(3*(i-1)+1, :);
        base = median(rr(i, tr <= 30));
        k = find(Ge > 3.6, 1);
        if isempty(k), g36 = NaN; else, g36 = k/cfg.eeg_fs; end
        above = rr(i,:) > max(4*base, 2.0);
        te = NaN;
        for j = 1:nr-5
            if all(above(j:j+4)), te = tr(j)-0.5; break; end
        end
        ref = make_reference(y(i, 1:round(10*fs)), fs, p);
        [~, labs] = classify_timeseries_v2(y(i,:), fs, p, ref);
        fo = sum(strcmp(labs,'ONSET'))/numel(labs);
        fi = sum(strcmp(labs,'ICTAL'))/numel(labs);
        fprintf('%-7.2f %-5s %9.3f %8.2f %8.2f %9.1f %8.1f %10.2f %10.2f\n', ...
            wlist(wi), nm{i}, base, max(rr(i,:)), max(Ge), g36, te, fo, fi);
        fprintf(fid, '%-7.2f %-5s %9.3f %8.2f %8.2f %9.1f %8.1f %10.2f %10.2f\n', ...
            wlist(wi), nm{i}, base, max(rr(i,:)), max(Ge), g36, te, fo, fi);
    end
end
fclose(fid);
fprintf('\n已写 results/niz_weight.txt\n');
end
